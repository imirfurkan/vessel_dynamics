// Dynamics implementation of the Milliampere vessel from NTNU

#include <memory>
#include <cmath>
#include <rclcpp/rclcpp.hpp>
#include "vessel_kinematics/thrust_allocator.hpp"

#include <geometry_msgs/msg/wrench.hpp>
#include <nav_msgs/msg/odometry.hpp>
#include <geometry_msgs/msg/transform_stamped.hpp>
#include <tf2_ros/transform_broadcaster.h>
#include <tf2/LinearMath/Quaternion.h>
#include <std_msgs/msg/float64_multi_array.hpp>

#include <Eigen/Dense>

using namespace std::chrono_literals;

class MilliampereDynamics : public rclcpp::Node
{
public:
  MilliampereDynamics()              // constructor of the class
      : Node("milliampere_dynamics") // initializer list that calls the base class class constructor
                                     // (rclcpp::Node) with the name "milliampere_dynamics"/
  {
    //----------------------------------------------------------------
    // 1) Parameters
    //----------------------------------------------------------------
    // clang-format off
    M_ << m11_,    0,    0,
             0, m22_, m23_,
             0, m32_, m33_;
    // clang-format on

    //----------------------------------------------------------------
    // 2) State initialization
    //----------------------------------------------------------------
    nu_.setZero();  // [u, v, r]
    eta_.setZero(); // [x, y, ψ]
    tau_.setZero(); // [X, Y, N]
    tau_des_.setZero();
    azimuthAngles_.setZero(); // radian
    // azimuthAngles_ << 10.0, 10.0;
    thrustForces_.setZero();

    //----------------------------------------------------------------
    // 3) ROS interfaces
    //----------------------------------------------------------------
    tau_cmd_sub_ = create_subscription<geometry_msgs::msg::Wrench>(
        "/tau_cmd",
        10,
        std::bind(&MilliampereDynamics::onTauCommand, this, std::placeholders::_1)); // TODO anla

    odom_pub_ = create_publisher<nav_msgs::msg::Odometry>("/odom", 10);

    // coordinate frame transforms (map → base_link)
    tf_broadcaster_ = std::make_shared<tf2_ros::TransformBroadcaster>(this);

    // calls step() each 20ms
    timer_ = create_wall_timer(20ms, std::bind(&MilliampereDynamics::step, this));
  }

private:
  // callback: capture incoming tau
  void onTauCommand(const geometry_msgs::msg::Wrench::SharedPtr msg)
  {
    tau_des_(0) = msg->force.x;  // desired X force
    tau_des_(1) = msg->force.y;  // desired Y force
    tau_des_(2) = msg->torque.z; // desired N torque
  }

  void step()
  {
    const double dt = 0.02; // 50 Hz

    // ------------------------------------------------------------
    //  Thruster Allocation
    // ------------------------------------------------------------
    RCLCPP_INFO(
        this->get_logger(), "Desired Tau: [%f, %f, %f]", tau_des_(0), tau_des_(1), tau_des_(2));
    TAResult ta_result = allocate_tau(tau_des_, ta_state_, ta_params_, dt);

    if (ta_result.success)
    {
      azimuthAngles_ = ta_result.alpha;
      thrustForces_  = ta_result.f;
    }
    else
    {
      // If the allocator fails, command zero thrust as a safety measure.
      thrustForces_.setZero();
      RCLCPP_WARN(this->get_logger(), "Thrust allocator failed to find a solution!");
    }

    ta_state_.alpha = ta_result.alpha;
    ta_state_.f     = ta_result.f;

    // ------------------------------------------------------------
    //  Rigid‐body dynamics
    // ------------------------------------------------------------
    double ca1 = std::cos(azimuthAngles_(0)), ca2 = std::cos(azimuthAngles_(1));
    double sa1 = std::sin(azimuthAngles_(0)), sa2 = std::sin(azimuthAngles_(1));

    // clang-format off
    T_ <<       ca1,      ca2,
                sa1,      sa2,
          Lx_ * sa1, -Lx_ * sa2; //
    // clang-format on

    T_   = Tmatrix(azimuthAngles_, Lx_);
    tau_ = T_ * thrustForces_;

    // debug
    RCLCPP_INFO(this->get_logger(),
                "Azimuths: [%.2f, %.2f], Thrusts: [%.2f, %.2f]",
                azimuthAngles_(0),
                azimuthAngles_(1),
                thrustForces_(0),
                thrustForces_(1));
    RCLCPP_INFO(
        this->get_logger(), "Actual Tau: [X=%.2f, Y=%.2f, N=%.2f]", tau_(0), tau_(1), tau_(2));

    // Mν̇ + C(v)v + D(v)v = τ
    updateCoriolisMatrix(nu_);
    updateDampingMatrix(nu_);
    Eigen::Vector3d nu_dot = M_.ldlt().solve(-C_ * nu_ - D_ * nu_ + tau_);
    // RCLCPP_INFO(get_logger(), "C*nu: [%f, %f, %f]", (C_ * nu_)(0), (C_ * nu_)(1), (C_ * nu_)(2));
    // RCLCPP_INFO(get_logger(), "D*nu: [%f, %f, %f]", (D_ * nu_)(0), (D_ * nu_)(1), (D_ * nu_)(2));
    // RCLCPP_INFO(get_logger(), "tau:  [%f, %f, %f]", tau_(0), tau_(1), tau_(2));

    // integrate body‐frame velocity νk+1​=νk​+ν˙k​⋅Δt
    nu_ += nu_dot * dt;

    RCLCPP_INFO(this->get_logger(), "nu: [u=%.2f, v=%.2f, r=%.2f]", nu_(0), nu_(1), nu_(2));
    // RCLCPP_INFO(
    //     this->get_logger(), "nu_dot: [du=%.3f, dv=%.3f, dr=%.3f]", nu_dot(0), nu_dot(1),
    //     nu_dot(2));

    // ------------------------------------------------------------
    //  Kinematic mapping to inertial pose
    // ------------------------------------------------------------

    double cy = std::cos(eta_(2)), sy = std::sin(eta_(2));

    // clang-format off
    Jnb << cy, -sy, 0,
           sy,  cy, 0,
            0,   0, 1;
    // clang-format on

    Eigen::Vector3d eta_dot = Jnb * nu_;
    eta_ += eta_dot * dt;

    // ------------------------------------------------------------
    //  Publish TF + Odometry
    // ------------------------------------------------------------
    publishState(this->now());
  }

  void publishState(const rclcpp::Time& stamp) // cpp note: for member functions inside a class, the
                                               // order of their definitions does not matter.
  {
    // 1) TF: map → base_link
    geometry_msgs::msg::TransformStamped tf_msg;
    tf_msg.header.stamp    = stamp;
    tf_msg.header.frame_id = "map";
    tf_msg.child_frame_id  = "base_link";
    // ---- compute ENU pose from NED state (2D) ----
    const double x_ros   = eta_(1);              // ENU x = East  = y_NED
    const double y_ros   = eta_(0);              // ENU y = North = x_NED
    const double yaw_ros = M_PI / 2.0 - eta_(2); // ENU yaw CCW = 90° - NED yaw (CW+)

    tf_msg.transform.translation.x = x_ros;
    tf_msg.transform.translation.y = y_ros;
    tf2::Quaternion q;
    q.setRPY(0.0, 0.0, yaw_ros);
    tf_msg.transform.rotation.x = q.x();
    tf_msg.transform.rotation.y = q.y();
    tf_msg.transform.rotation.z = q.z();
    tf_msg.transform.rotation.w = q.w();
    tf_broadcaster_->sendTransform(tf_msg);

    // 2) Odometry
    nav_msgs::msg::Odometry odom;
    odom.header.stamp    = stamp;
    odom.header.frame_id = "map";
    odom.child_frame_id  = "base_link";

    // pose
    odom.pose.pose.position.x  = x_ros;
    odom.pose.pose.position.y  = y_ros;
    odom.pose.pose.position.z  = 0.0; // flat world
    odom.pose.pose.orientation = tf_msg.transform.rotation;

    // twist (body frame)
    odom.twist.twist.linear.x  = nu_(0);
    odom.twist.twist.linear.y  = nu_(1);                //
    odom.twist.twist.angular.z = nu_(2) * 180.0 / M_PI; // revert to normal later.

    odom_pub_->publish(odom);
  }

  // ---------- ROS interfaces --------------------- // TODO anla
  rclcpp::Subscription<geometry_msgs::msg::Wrench>::SharedPtr tau_cmd_sub_;
  rclcpp::Publisher<nav_msgs::msg::Odometry>::SharedPtr       odom_pub_;
  std::shared_ptr<tf2_ros::TransformBroadcaster>              tf_broadcaster_;
  rclcpp::TimerBase::SharedPtr                                timer_;

  // ---------- node state & parameters ------------
  double          m_  = 1800;     // mass [kg]
  double          Lx_ = 1.8;      // distance from CO to thrusters [m] TODO find the actual value
  Eigen::Vector3d tau_des_;       // [X, Y, N]
  Eigen::Matrix3d Jnb;            // rotation matrix from body to inertial
  Eigen::Matrix3d M_;             // mass inertia matrix
  Eigen::Matrix3d C_;             // coriolis-centripetal matrix
  Eigen::Matrix3d D_;             // damping matrix
  Eigen::Vector3d eta_;           // position [x, y, ψ]
  Eigen::Vector3d nu_;            // body-frame twist [u, v, r]
  Eigen::Vector3d tau_;           // [X, Y, N]
  Eigen::Vector2d azimuthAngles_; // [α₁, α₂] in radians
  Eigen::Vector2d thrustForces_;  // [F₁, F₂] in Newtons
  Eigen::Matrix<double, 3, 2> T_; // Thrust configuration matrix

  TAState  ta_state_;
  TAParams ta_params_;

  double m11_ = 2389.657; // [kg]
  double m22_ = 2533.911; // [kg]
  double m23_ = 62.386;   // [kg]
  double m32_ = 28.141;   // [kg]
  double m33_ = 5068.910; // [kg] or [kg*m2]

  // ---------- Damping matrix D(nu) coefficients (fully coupled model) ----------
  // Surge
  double Xu_   = -27.632;  // [kg/s]
  double Xuu_  = -110.064; // [kg/s]
  double Xuuu_ = -13.965;  // [kg/s]
  // Sway
  double Yv_   = -52.947;   // [kg/s]
  double Yvv_  = -116.486;  // [kg/s]
  double Yvvv_ = -24.313;   // [kg/s]
  double Yvr_  = 572.141;   // [kg/s]
  double Yrv_  = -1540.383; // [kg/s]
  double Yrr_  = -115.457;  // [kg/s]
  double Yr_   = 24.372;    // [kg/s]
  // Yaw
  double Nv_   = 3.524;    // [kg/s]
  double Nrv_  = 336.827;  // [kg/s]
  double Nvv_  = -0.832;   // [kg/s]
  double Nr_   = -122.860; // [kg/s]
  double Nrr_  = -874.428; // [kg/s]
  double Nvr_  = -121.957; // [kg/s]
  double Nrrr_ = 0.0;      // [kg/s]

  // velocity dependent → compute dynamically
  void updateCoriolisMatrix(const Eigen::Vector3d& nu)
  {
    double u = nu(0), v = nu(1), r = nu(2);

    double c13 = -m22_ * v - m23_ * r;
    double c23 = m11_ * u;

    // clang-format off
    C_ <<  0.0,  0.0, c13,
           0.0,  0.0, c23,
          -c13, -c23, 0.0;
    // clang-format on
  }

  // velocity dependent → compute dynamically
  void updateDampingMatrix(const Eigen::Vector3d& nu)
  {
    double u = nu(0), v = nu(1), r = nu(2);

    double d11 = -Xu_ - Xuu_ * std::abs(u) - Xuuu_ * u * u;
    double d22 = -Yv_ - Yvv_ * std::abs(v) - Yrv_ * std::abs(r) - Yvvv_ * v * v;
    double d23 = -Yr_ - Yvr_ * std::abs(v) - Yrr_ * std::abs(r);
    double d32 = -Nv_ - Nvv_ * std::abs(v) - Nrv_ * std::abs(r);
    double d33 = -Nr_ - Nvr_ * std::abs(v) - Nrr_ * std::abs(r) - Nrrr_ * r * r;

    // clang-format off
    D_ << d11, 0.0,  0.0,
          0.0, d22,  d23,
          0.0, d32,  d33;
    // clang-format on
  }

  double deg2rad(double deg) { return deg * M_PI / 180.0; }

  double wrapAngle(double angle)
  {
    const double two_pi = 2.0 * M_PI;

    // Step 1: First, wrap the angle into the [0, 2π) range.
    // After this block, 'angle' is guaranteed to be positive.
    angle = std::fmod(angle, two_pi);
    if (angle < 0)
      angle += two_pi;

    // Step 2: Now, move angles from the top half (π, 2π) to the equivalent
    // negative half (-π, 0).
    if (angle > M_PI)
      angle -= two_pi;

    return angle;
  }
};

int main(int argc, char** argv)
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<MilliampereDynamics>());
  rclcpp::shutdown();
  return 0;
}