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

using namespace std::chrono_literals; // write time as ms

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
    tau_wave_.setZero();

    //----------------------------------------------------------------
    // ROS interfaces
    //----------------------------------------------------------------
    actual_wrench_sub_ = create_subscription<geometry_msgs::msg::Wrench>(
        "/actual_wrench", 10, std::bind(&MilliampereDynamics::actualWrenchCallback, this, std::placeholders::_1));

    disturbance_sub_ = create_subscription<geometry_msgs::msg::Wrench>(
        "/wave_wrench",
        10,
        std::bind(&MilliampereDynamics::waveWrenchCallback, this, std::placeholders::_1)); // TODO anla

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

  void waveWrenchCallback(const geometry_msgs::msg::Wrench::SharedPtr msg)
  {
    tau_wave_(0) = msg->force.x;
    tau_wave_(1) = msg->force.y;
    tau_wave_(2) = msg->torque.z;
  }

  void step()
  {
    // ------------------------------------------------------------
    //  Rigid‐body dynamics -> Mν̇ + C(v)v + D(v)v = τ + t_wave
    // ------------------------------------------------------------

    Eigen::Vector3d tau_total_ = tau_actual_ + tau_wave_;

    vk::updateCoriolisMatrix(nu_, this->C_); // TODO understand this->
    vk::updateDampingMatrix(nu_, this->D_);
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
    RCLCPP_INFO(this->get_logger(), "Position: [N=%.2f, E=%.2f, ψ=%.2f°]", eta_(0), eta_(1), eta_(2) * 180 / M_PI);
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
  rclcpp::Subscription<geometry_msgs::msg::Wrench>::SharedPtr      actual_wrench_sub_;
  rclcpp::Subscription<geometry_msgs::msg::Wrench>::SharedPtr      disturbance_sub_;
  rclcpp::Publisher<nav_msgs::msg::Odometry>::SharedPtr            odom_pub_;
  std::shared_ptr<tf2_ros::TransformBroadcaster>                   tf_broadcaster_;
  rclcpp::TimerBase::SharedPtr                                     timer_;

  // ---------- node state & parameters ------------
  const double    m_  = 1800;  // mass [kg]
  const double    Lx_ = 1.8;   // distance from CO to thrusters [m] TODO find the actual value
  Eigen::Vector3d eta_;        // position [x, y, ψ]
  Eigen::Vector3d nu_;         // body-frame twist [u, v, r]
  Eigen::Vector3d tau_actual_; // [X, Y, N]
  Eigen::Vector3d tau_wave_;   // [X, Y, N]
  Eigen::Matrix3d Jnb;         // rotation matrix from body to inertial
  Eigen::Matrix3d C_;          // coriolis-centripetal matrix
  Eigen::Matrix3d D_;          // damping matrix
};

int main(int argc, char** argv)
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<MilliampereDynamics>());
  rclcpp::shutdown();
  return 0;
}