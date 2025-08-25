// Dynamics implementation of the Milliampere vessel from NTNU

#include <memory>
#include <cmath>
#include <rclcpp/rclcpp.hpp>
#include "vessel_kinematics/thrust_allocator.hpp"
#include <geometry_msgs/msg/wrench.hpp>
#include <nav_msgs/msg/odometry.hpp>
#include <geometry_msgs/msg/transform_stamped.hpp>
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <tf2_ros/transform_broadcaster.h>
#include <tf2/LinearMath/Quaternion.h>
#include <std_msgs/msg/float64_multi_array.hpp>
#include <tf2/LinearMath/Matrix3x3.h>
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

    Kp.diagonal() << 200, 200, 800;
    Ki.diagonal() << 10, 10, 45;
    Kd.diagonal() << 700, 700, 1600;

    integral_error_limits_ << 50.0, 50.0, 50.0; // TODO tune these

    //----------------------------------------------------------------
    // 2) State initialization
    //----------------------------------------------------------------
    nu_.setZero();  // [u, v, r]
    eta_.setZero(); // [x, y, ψ]
    tau_.setZero(); // [X, Y, N]
    tau_des_.setZero();
    nu_des_.setZero();
    azimuthAngles_.setZero(); // radian
    thrustForces_.setZero();
    eta_ref_ = eta_; // Start the reference model at the vessel's initial state
    eta_dot_ref_.setZero();
    eta_ddot_ref_.setZero();

    vel_ref_max_ << 2.57,    // Max speed ~ 5 knots [m/s]
        1.0,                 // TODO (arbitrary)
        32.0 * M_PI / 180.0; // Table 5.1 (not azimuth but body?)

    acc_ref_max_ << 0.25, // TODO (arbitrary)
        0.1,              // TODO (arbitrary)
        0.2;              // TODO (arbitrary)

    // tuning parameters for the reference model
    // reference model should be a slow, smooth system that ignores fast disturbances
    // as we want reference feedforward model to be low bandwidth and pid to be high bandwidth
    ref_zeta_ << 1.0, 1.0, 1.0;     // TODO Critically damped is usually best
    ref_omega_n_ << 0.5, 0.25, 0.5; // TODO Tune these for response speed, higher values create
                                    // high-bandwidth reference signal

    //----------------------------------------------------------------
    // 3) ROS interfaces
    //----------------------------------------------------------------
    pos_cmd_sub_ = create_subscription<geometry_msgs::msg::PoseStamped>(
        "/goal_pose", 10, std::bind(&MilliampereDynamics::poseCallback, this, std::placeholders::_1)); // TODO anla

    odom_pub_ = create_publisher<nav_msgs::msg::Odometry>("/odom", 10);

    // coordinate frame transforms (map → base_link)
    tf_broadcaster_ = std::make_shared<tf2_ros::TransformBroadcaster>(this);

    // calls step() each 20ms
    timer_       = create_wall_timer(20ms, std::bind(&MilliampereDynamics::step, this));
    azimuth_pub_ = create_publisher<std_msgs::msg::Float64MultiArray>("/azimuth_angles", 10);
  }

private:
  // callback: capture incoming tau
  void poseCallback(const geometry_msgs::msg::PoseStamped::SharedPtr msg)
  {
    // tf2::Quaternion Q(
    //     msg->pose.orientation.x, msg->pose.orientation.y, msg->pose.orientation.z, msg->pose.orientation.w);
    // tf2::Matrix3x3 m(Q);
    // double         roll, pitch, yaw;
    // m.getRPY(roll, pitch, yaw);

    double yaw_deg = msg->pose.orientation.z;
    double yaw_rad = yaw_deg * M_PI / 180;

    // Be careful when converting from ROS convention ENU to NED frames.
    eta_des_(0) = msg->pose.position.y;
    eta_des_(1) = msg->pose.position.x;
    eta_des_(2) = yaw_rad;

    RCLCPP_INFO(this->get_logger(), "Reference Position: [%f, %f, %f]", eta_des_(0), eta_des_(1), eta_des_(2));
  }

  void step()
  {
    const double dt = 0.02; // 50 Hz

    // ------------------------------------------------------------
    //  Reference Model + FF + PID
    // ------------------------------------------------------------
    Eigen::Matrix3d Rbn = calculateRotationMatrix(eta_).transpose(); // of eta_, not eta_des_

    updateReferenceModel(dt);
    Eigen::Vector3d nu_ref_     = Rbn * eta_dot_ref_;  // velocity in the body frame
    Eigen::Vector3d nu_dot_ref_ = Rbn * eta_ddot_ref_; // acceleration in the body frame
    RCLCPP_INFO(this->get_logger(), "nu_ref_: [%f, %f, %f]", nu_ref_(0), nu_ref_(1), nu_ref_(2));
    RCLCPP_INFO(this->get_logger(), "nu_dot_ref_: [%f, %f, %f]", nu_dot_ref_(0), nu_dot_ref_(1), nu_dot_ref_(2));

    Eigen::Matrix3d C_d, D_d;
    updateCoriolisMatrix(nu_ref_, C_d);
    updateDampingMatrix(nu_ref_, D_d);
    Eigen::Vector3d tau_ff = M_ * nu_dot_ref_ + C_d * nu_ref_ + D_d * nu_ref_;

    Eigen::Vector3d eta_diff_;
    eta_diff_ << (eta_(0) - eta_ref_(0)), // N error
        (eta_(1) - eta_ref_(1)),          // E error
        wrapAngle(eta_(2) - eta_ref_(2)); // wrapped ψ error

    Eigen::Vector3d error_     = Rbn * eta_diff_; // error should be in the body frame
    Eigen::Vector3d error_dot_ = nu_ - nu_ref_;

    error_integral_ += error_ * dt;
    error_integral_ = clampVec3(error_integral_, -integral_error_limits_, integral_error_limits_);

    Eigen::Vector3d tau_pid = -Kp * error_ - Ki * error_integral_ - Kd * error_dot_; // TODO

    tau_des_ = tau_pid + tau_ff; // TODO add disturbances
    RCLCPP_INFO(this->get_logger(), "FF Tau: [%f, %f, %f]", tau_ff(0), tau_ff(1), tau_ff(2));
    RCLCPP_INFO(this->get_logger(), "PID Tau: [%f, %f, %f]", tau_pid(0), tau_pid(1), tau_pid(2));
    RCLCPP_INFO(this->get_logger(), "Desired Tau: [%f, %f, %f]", tau_des_(0), tau_des_(1), tau_des_(2));

    // ------------------------------------------------------------
    //  Thruster Allocation
    // ------------------------------------------------------------
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

    // Create the message object
    auto msg = std::make_shared<std_msgs::msg::Float64MultiArray>();

    // Resize the data to hold 2 elements
    msg->data.resize(4);

    // Copy the Eigen::Vector2d data into the message data vector
    msg->data[0] = azimuthAngles_(0) * 180 / M_PI;
    msg->data[1] = azimuthAngles_(1) * 180 / M_PI;
    msg->data[2] = tau_pid(0);
    msg->data[3] = tau_ff(0);

    // Publish the message
    azimuth_pub_->publish(*msg);

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
                "Azimuths: [%.2f, %.2f], Thrusts: [%.2f, %.2f]", // TODO change output to degrees
                azimuthAngles_(0) * 180 / M_PI,
                azimuthAngles_(1) * 180 / M_PI,
                thrustForces_(0),
                thrustForces_(1));
    RCLCPP_INFO(this->get_logger(), "Actual Tau: [X=%.2f, Y=%.2f, N=%.2f]", tau_(0), tau_(1), tau_(2));
    RCLCPP_INFO(this->get_logger(), "Position: [N=%.2f, E=%.2f, ψ=%.2f]", eta_(0), eta_(1), eta_(2));

    // Mν̇ + C(v)v + D(v)v = τ

    updateCoriolisMatrix(nu_, this->C_);
    updateDampingMatrix(nu_, this->D_);
    Eigen::Vector3d nu_dot = M_.ldlt().solve(-C_ * nu_ - D_ * nu_ + tau_);
    // RCLCPP_INFO(get_logger(), "C*nu: [%f, %f, %f]", (C_ * nu_)(0), (C_ * nu_)(1), (C_ * nu_)(2));
    // RCLCPP_INFO(get_logger(), "D*nu: [%f, %f, %f]", (D_ * nu_)(0), (D_ * nu_)(1), (D_ * nu_)(2));
    // RCLCPP_INFO(get_logger(), "tau:  [%f, %f, %f]", tau_(0), tau_(1), tu_(2));a

    nu_ += nu_dot * dt;

    RCLCPP_INFO(this->get_logger(), "nu: [u=%.2f, v=%.2f, r=%.2f]", nu_(0), nu_(1), nu_(2));

    // ------------------------------------------------------------
    //  Kinematic mapping to inertial pose
    // ------------------------------------------------------------

    Eigen::Matrix3d Rnb     = calculateRotationMatrix(eta_);
    Eigen::Vector3d eta_dot = Rnb * nu_;
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
  rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr pos_cmd_sub_;
  rclcpp::Publisher<nav_msgs::msg::Odometry>::SharedPtr            odom_pub_;
  std::shared_ptr<tf2_ros::TransformBroadcaster>                   tf_broadcaster_;
  rclcpp::TimerBase::SharedPtr                                     timer_;
  rclcpp::Publisher<std_msgs::msg::Float64MultiArray>::SharedPtr   azimuth_pub_;

  // ---------- node state & parameters ------------
  const double                m_  = 1800;        // mass [kg]
  const double                Lx_ = 1.8;         // distance from CO to thrusters [m] TODO find the actual value
  Eigen::Vector3d             eta_;              // position [x, y, ψ]
  Eigen::Vector3d             nu_;               // body-frame twist [u, v, r]
  Eigen::Vector3d             tau_;              // [X, Y, N]
  Eigen::Vector3d             eta_des_;          // [N, E, ψ]
  Eigen::Vector3d             nu_des_;           // desired body-frame twist [u, v, r]
  Eigen::Vector3d             tau_des_;          // [X, Y, N]
  Eigen::Vector3d             eta_ref_;          // η_d   (smooth desired position)
  Eigen::Vector3d             eta_dot_ref_;      // η̇_d or v_d (smooth desired velocity)
  Eigen::Vector3d             eta_ddot_ref_;     // η̈_d or a_d (smooth desired acceleration)
  Eigen::Matrix3d             Jnb;               // rotation matrix from body to inertial
  Eigen::Matrix3d             M_;                // mass inertia matrix
  Eigen::Matrix3d             C_;                // coriolis-centripetal matrix
  Eigen::Matrix3d             D_;                // damping matrix
  Eigen::Matrix3d             Kp;                //
  Eigen::Matrix3d             Ki;                //
  Eigen::Matrix3d             Kd;                //
  Eigen::Vector3d             error_integral_{}; //
  Eigen::Vector3d             integral_error_limits_;
  Eigen::Vector2d             azimuthAngles_; // [α₁, α₂] in radians
  Eigen::Vector2d             thrustForces_;  // [F₁, F₂] in Newtons
  Eigen::Matrix<double, 3, 2> T_;             // Thrust configuration matrix

  // Physical limits for saturation
  Eigen::Vector3d vel_ref_max_; // v_max in the book
  Eigen::Vector3d acc_ref_max_; // a_max (derived from ν̇_max)

  // Tuning parameters for the reference model (Ω and Δ matrices)
  Eigen::Vector3d ref_omega_n_; // Diagonal elements of Ω
  Eigen::Vector3d ref_zeta_;    // Diagonal elements of Δ

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
  void updateCoriolisMatrix(const Eigen::Vector3d& nu, Eigen::Matrix3d& C_out)
  {
    double u = nu(0), v = nu(1), r = nu(2);

    double c13 = -m22_ * v - m23_ * r;
    double c23 = m11_ * u;

    // clang-format off
      C_out <<  0.0,  0.0, c13,
            0.0,  0.0, c23,
            -c13, -c23, 0.0;
    // clang-format on
  }

  // velocity dependent → compute dynamically
  void updateDampingMatrix(const Eigen::Vector3d& nu, Eigen::Matrix3d& D_out)
  {
    double u = nu(0), v = nu(1), r = nu(2);

    double d11 = -Xu_ - Xuu_ * std::abs(u) - Xuuu_ * u * u;
    double d22 = -Yv_ - Yvv_ * std::abs(v) - Yrv_ * std::abs(r) - Yvvv_ * v * v;
    double d23 = -Yr_ - Yvr_ * std::abs(v) - Yrr_ * std::abs(r);
    double d32 = -Nv_ - Nvv_ * std::abs(v) - Nrv_ * std::abs(r);
    double d33 = -Nr_ - Nvr_ * std::abs(v) - Nrr_ * std::abs(r) - Nrrr_ * r * r;

    // clang-format off
      D_out << d11, 0.0,  0.0,
            0.0, d22,  d23,
            0.0, d32,  d33;
    // clang-format on
  }

  Eigen::Matrix3d calculateRotationMatrix(const Eigen::Vector3d& pos)
  {
    double cy = std::cos(pos(2)), sy = std::sin(pos(2)); // pos(2) is yaw.

    // clang-format off
      Jnb << cy, -sy, 0,
            sy,  cy, 0,
              0,   0, 1;
    // clang-format on

    return Jnb;
  }

  void updateReferenceModel(double dt)
  {
    // Let w = ref_omega_n_ and z = ref_zeta_ for simplicity.
    // The equation is: η3_d + (2Δ+I)Ωη̈_d + (2Δ+I)Ω²η̇_d + Ω³η_d = Ω³rⁿ
    // We solve for the highest derivative, η3_d (desired jerk).

    // Diagonals of Ω matrices
    Eigen::Vector3d w  = ref_omega_n_;                            // Ω
    Eigen::Vector3d w2 = ref_omega_n_.cwiseProduct(ref_omega_n_); // Ω²
    Eigen::Vector3d w3 = w2.cwiseProduct(ref_omega_n_);           // Ω³

    // Diagonals of the (2Δ + I) matrix
    Eigen::Vector3d two_z_plus_I = 2.0 * ref_zeta_ + Eigen::Vector3d::Ones();

    // Calculate desired jerk by building the rearranged equation term by term
    Eigen::Vector3d term1 =
        w3.cwiseProduct(eta_des_ - eta_ref_); // TODO does it make sense for eta_ref_ to be initialized with eta_
    Eigen::Vector3d term2 = w2.cwiseProduct(two_z_plus_I).cwiseProduct(eta_dot_ref_);
    Eigen::Vector3d term3 = w.cwiseProduct(two_z_plus_I).cwiseProduct(eta_ddot_ref_);

    Eigen::Vector3d eta_dddot_ref_ = term1 - term2 - term3;

    eta_ddot_ref_ += eta_dddot_ref_ * dt; // Eular integration from jerk to acceleration
    eta_ddot_ref_ = clampVec3(eta_ddot_ref_, -acc_ref_max_, acc_ref_max_); // Saturate the acceleration
    eta_dot_ref_ += eta_ddot_ref_ * dt;
    eta_dot_ref_ = clampVec3(eta_dot_ref_, -vel_ref_max_, vel_ref_max_);
    eta_ref_ += eta_dot_ref_ * dt;
    eta_ref_(2) = wrapAngle(eta_ref_(2));
  }

  double deg2rad(double deg) { return deg * M_PI / 180.0; }

  double wrapAngle(double angle)
  {
    const double two_pi = 2.0 * M_PI;

    angle = std::fmod(angle, two_pi);
    if (angle < 0)
      angle += two_pi;

    if (angle > M_PI)
      angle -= two_pi;

    return angle;
  }

  Eigen::Vector3d clampVec3(const Eigen::Vector3d& x, const Eigen::Vector3d& lo, const Eigen::Vector3d& hi)
  {
    return x.cwiseMax(lo).cwiseMin(hi);
  }
};

int main(int argc, char** argv)
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<MilliampereDynamics>());
  rclcpp::shutdown();
  return 0;
}