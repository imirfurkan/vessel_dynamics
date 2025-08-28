// Dynamics implementation of the Milliampere vessel from NTNU

#include <memory>
#include <cmath>
#include <rclcpp/rclcpp.hpp>
#include "vessel_kinematics/thrust_allocator.hpp"
#include "vessel_kinematics/utils.hpp"
#include <geometry_msgs/msg/wrench.hpp>
#include <nav_msgs/msg/odometry.hpp>
#include <geometry_msgs/msg/transform_stamped.hpp>
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <tf2_ros/transform_broadcaster.h>
#include <tf2/LinearMath/Quaternion.h>
#include <std_msgs/msg/float64_multi_array.hpp>
#include <tf2/LinearMath/Matrix3x3.h>
#include <Eigen/Dense>
#include <random>
#include <chrono>

using namespace std::chrono_literals;

class MilliampereDynamics : public rclcpp::Node
{
public:
  MilliampereDynamics()              // constructor of the class
      : Node("milliampere_dynamics") // initializer list that calls the base class class constructor
                                     // (rclcpp::Node) with the name "milliampere_dynamics"/
  {
    //----------------------------------------------------------------
    // State initialization
    //----------------------------------------------------------------
    nu_.setZero();         // [u, v, r]
    eta_.setZero();        // [x, y, ψ]
    tau_actual_.setZero(); // [X, Y, N]

    // Wave states
    omega0_ << 0.90, 0.90, 0.90; // rad/s
    lambda_ << 0.10, 0.10, 0.10; // JONSWAP damping ratio // TODO why jonswap
    K_wave_ << 500.0, 500.0, 500.0;
    xi_wave_.setZero();

    A_wave_.setZero();
    A_wave_.topRightCorner<3, 3>() = Eigen::Matrix3d::Identity();

    Eigen::Matrix3d omega0_sq        = omega0_.cwiseProduct(omega0_).asDiagonal();
    A_wave_.bottomLeftCorner<3, 3>() = -omega0_sq;

    Eigen::Matrix3d TwoLambdaOmega0   = (2.0 * lambda_.array() * omega0_.array()).matrix().asDiagonal();
    A_wave_.bottomRightCorner<3, 3>() = -TwoLambdaOmega0;

    E_wave_.setZero();
    E_wave_.bottomRows<3>() = Eigen::Matrix3d::Identity();

    //----------------------------------------------------------------
    // 3) ROS interfaces
    //----------------------------------------------------------------
    actual_wrench_sub_ = create_subscription<geometry_msgs::msg::Wrench>(
        "/actual_wrench", 10, std::bind(&MilliampereDynamics::actualWrenchCallback, this, std::placeholders::_1));

    odom_pub_ = create_publisher<nav_msgs::msg::Odometry>("/vessel_state", 10);

    // coordinate frame transforms (map → base_link)
    tf_broadcaster_ = std::make_shared<tf2_ros::TransformBroadcaster>(this);

    // calls step() each 20ms
    timer_ = create_wall_timer(20ms, std::bind(&MilliampereDynamics::step, this));
  }

private:
  // callback: capture incoming tau
  void actualWrenchCallback(const geometry_msgs::msg::Wrench::SharedPtr msg)
  {
    tau_actual_(0) = msg->force.x;
    tau_actual_(1) = msg->force.y;
    tau_actual_(2) = msg->torque.z;
  }

  void step()
  {
    // TODO add msg from the other node

    // ------------------------------------------------------------
    //  Rigid‐body dynamics
    // ------------------------------------------------------------

    // Mν̇ + C(v)v + D(v)v = τ + t_wave

    std::default_random_engine       generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::normal_distribution<double> distribution(0.0, 1.0); // Mean 0, StdDev 1

    Eigen::Vector3d white_noise_ = {distribution(generator), distribution(generator), distribution(generator)};

    Eigen::Vector<double, 6> xi_dot = A_wave_ * xi_wave_ + E_wave_ * white_noise_;
    xi_wave_ += xi_dot * vk::dt;

    Eigen::Vector3d first_order_waves = K_wave_.cwiseProduct(xi_wave_.tail<3>());
    tau_wave_                         = first_order_waves; // TODO + second_order_waves

    vk::updateCoriolisMatrix(nu_, this->C_); // TODO understand this->
    vk::updateDampingMatrix(nu_, this->D_);
    Eigen::Vector3d tau_total_ = tau_actual_ + tau_wave_;
    // Eigen::Vector3d tau_total_ = tau_wave_;

    Eigen::Vector3d nu_dot = vk::M_.ldlt().solve(-C_ * nu_ - D_ * nu_ + tau_total_);
    // RCLCPP_INFO(get_logger(), "C*nu: [%f, %f, %f]", (C_ * nu_)(0), (C_ * nu_)(1), (C_ * nu_)(2));
    // RCLCPP_INFO(get_logger(), "D*nu: [%f, %f, %f]", (D_ * nu_)(0), (D_ * nu_)(1), (D_ * nu_)(2));
    // RCLCPP_INFO(get_logger(), "tau:  [%f, %f, %f]", tau_(0), tau_(1), tu_(2));a

    nu_ += nu_dot * vk::dt;

    RCLCPP_INFO(this->get_logger(), "nu: [u=%.2f, v=%.2f, r=%.2f]", nu_(0), nu_(1), nu_(2));

    // ------------------------------------------------------------
    //  Kinematic mapping to inertial pose
    // ------------------------------------------------------------

    Eigen::Matrix3d Rnb     = vk::calculateRotationMatrix(eta_);
    Eigen::Vector3d eta_dot = Rnb * nu_;
    eta_ += eta_dot * vk::dt;
    RCLCPP_INFO(this->get_logger(), "Position: [N=%.2f, E=%.2f, ψ=%.2f° ]", eta_(0), eta_(1), eta_(2));
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
    const double x_ros   = eta_(1);                             // ENU x = East  = y_NED
    const double y_ros   = eta_(0);                             // ENU y = North = x_NED
    const double yaw_ros = vk::wrapAngle(M_PI / 2.0 - eta_(2)); // ENU yaw CCW = 90° - NED yaw (CW+)

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
    odom.twist.twist.linear.y  = nu_(1);
    odom.twist.twist.angular.z = nu_(2);

    odom_pub_->publish(odom);
  }

  // ---------- ROS interfaces ---------------------
  rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr pos_cmd_sub_;
  rclcpp::Publisher<nav_msgs::msg::Odometry>::SharedPtr            odom_pub_;
  std::shared_ptr<tf2_ros::TransformBroadcaster>                   tf_broadcaster_;
  rclcpp::TimerBase::SharedPtr                                     timer_;
  rclcpp::Subscription<geometry_msgs::msg::Wrench>::SharedPtr      actual_wrench_sub_;

  // ---------- node state & parameters ------------
  const double                m_  = 1800;  // mass [kg]
  const double                Lx_ = 1.8;   // distance from CO to thrusters [m] TODO find the actual value
  Eigen::Vector3d             eta_;        // position [x, y, ψ]
  Eigen::Vector3d             nu_;         // body-frame twist [u, v, r]
  Eigen::Vector3d             tau_actual_; // [X, Y, N]
  Eigen::Vector3d             tau_wave_;   // [X, Y, N]
  Eigen::Matrix3d             Jnb;         // rotation matrix from body to inertial
  Eigen::Matrix3d             C_;          // coriolis-centripetal matrix
  Eigen::Matrix3d             D_;          // damping matrix
  Eigen::Matrix<double, 3, 2> T_;          // Thrust configuration matrix

  TAState  ta_state_;
  TAParams ta_params_;

  // ---------- Wave Disturbance Model ------------
  Eigen::Vector3d             omega0_;  // Dominant wave frequencies [rad/s] for surge, sway, yaw
  Eigen::Vector3d             lambda_;  // Damping ratios for surge, sway, yaw
  Eigen::Vector3d             K_wave_;  // Gains to control force/torque magnitude
  Eigen::Vector<double, 6>    xi_wave_; // Grouped state vector [pos_s, pos_y, pos_n, vel_s, vel_y, vel_n]
  Eigen::Matrix<double, 6, 6> A_wave_;  // 6x6 state matrix
  Eigen::Matrix<double, 6, 3> E_wave_;  // 6x3 input matrix
};

int main(int argc, char** argv)
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<MilliampereDynamics>());
  rclcpp::shutdown();
  return 0;
}