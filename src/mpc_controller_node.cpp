#include <memory>
#include <rclcpp/rclcpp.hpp>
#include <nav_msgs/msg/odometry.hpp>
#include "geometry_msgs/msg/wrench.hpp"
#include <geometry_msgs/msg/pose_stamped.hpp>
#include "vessel_kinematics/utils.hpp"
#include <Eigen/Dense>
#include <chrono>
#include <casadi/casadi.hpp>
#include <tf2/LinearMath/Quaternion.h>
#include <tf2/LinearMath/Matrix3x3.h>

using namespace std::chrono_literals; // write time as ms
using namespace casadi;

class MPCNode : public rclcpp::Node
{

public:
  MPCNode() : Node("mpc_controller_node")
  {
    //----------------------------------------------------------------
    // ROS interfaces
    //----------------------------------------------------------------
    pos_cmd_sub_ = create_subscription<geometry_msgs::msg::PoseStamped>(
        "/goal_pose", 10, std::bind(&MPCNode::goalCallback, this, std::placeholders::_1));

    state_sub_ = create_subscription<nav_msgs::msg::Odometry>(
        "/vessel_state", 10, std::bind(&MPCNode::stateCallback, this, std::placeholders::_1));

    wrench_pub_ = create_publisher<geometry_msgs::msg::Wrench>("/cmd_wrench", 10);

    // calls step() each 20ms
    timer_ = create_wall_timer(20ms, std::bind(&MPCNode::controllerLoop, this));
  }

private:
  void goalCallback(const geometry_msgs::msg::PoseStamped::SharedPtr msg)
  {
    // NOTE: X and Y commands of goal_pose comply with ROS convention ENU frame,
    // while yaw is given directly as degrees in NED frame for convenience.
    goal_received_ = true;
    double yaw_deg = msg->pose.orientation.z; // command in rqt is given directly as degrees in NED frame.
    double yaw_rad = yaw_deg * M_PI / 180;

    // converting from ROS convention ENU to NED frames.
    eta_des_(0) = msg->pose.position.y;
    eta_des_(1) = msg->pose.position.x;
    eta_des_(2) = yaw_rad;

    RCLCPP_INFO(this->get_logger(), "Goal Position: [%f, %f, %f]", eta_des_(0), eta_des_(1), eta_des_(2));
  }

  void stateCallback(const nav_msgs::msg::Odometry::SharedPtr msg)
  {                                      // Odometry Twist
    nu_(0) = msg->twist.twist.linear.x;  // u
    nu_(1) = msg->twist.twist.linear.y;  // v
    nu_(2) = msg->twist.twist.angular.z; // r
    // Odometry Pose is published in ENU frame, convert to NED
    tf2::Quaternion q(msg->pose.pose.orientation.x,
                      msg->pose.pose.orientation.y,
                      msg->pose.pose.orientation.z,
                      msg->pose.pose.orientation.w);

    double roll, pitch, yaw_enu;
    tf2::Matrix3x3(q).getRPY(roll, pitch, yaw_enu);

    eta_(0) = msg->pose.pose.position.y;
    eta_(1) = msg->pose.pose.position.x;
    eta_(2) = vk::wrapAngle(M_PI / 2.0 - yaw_enu); //    // convert ENU yaw (CCW, 0 = East) to NED yaw (CW, 0 = North)
  }
  void controllerLoop() {}
  void setupSolver() {}
  // ROS2 Subscribers and Publishers
  rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr pos_cmd_sub_;
  rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr         state_sub_;
  rclcpp::Publisher<geometry_msgs::msg::Wrench>::SharedPtr         wrench_pub_;
  rclcpp::TimerBase::SharedPtr                                     timer_;

  // Variable declarations
  Eigen::Vector3d eta_;     // position [x, y, ψ]
  Eigen::Vector3d eta_des_; // [N, E, ψ]
  Eigen::Vector3d nu_;      // body-frame twist [u, v, r]

  bool goal_received_ = false;
};

int main(int argc, char** argv)
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<MPCNode>());
  rclcpp::shutdown();
  return 0;
}
