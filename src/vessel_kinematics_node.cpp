// src/vessel_kinematics_node.cpp

#include <cmath>
#include <memory>

#include <rclcpp/rclcpp.hpp>
#include <geometry_msgs/msg/twist.hpp>
#include <geometry_msgs/msg/transform_stamped.hpp>
#include <nav_msgs/msg/odometry.hpp>
#include <tf2_ros/transform_broadcaster.h>
#include <tf2/LinearMath/Quaternion.h>

using std::placeholders::_1;
using namespace std::chrono_literals;

class VesselKinematics : public rclcpp::Node
{
public:
  VesselKinematics()
  : Node("vessel_kinematics"),
    // initial pose
    x_(0.0), y_(0.0), z_(0.0),
    roll_(0.0), pitch_(0.0), yaw_(0.0),
    // initial body‐frame velocities
    u_(0.0), v_(0.0), w_(0.0),
    p_(0.0), q_(0.0), r_(0.0)
  {
    // subscriber to Twist in BODY frame
    sub_cmd_vel_ = create_subscription<geometry_msgs::msg::Twist>(
      "/cmd_vel_body", 10,
      std::bind(&VesselKinematics::cmdVelCallback, this, _1));

    // TF broadcaster
    tf_broadcaster_ = std::make_shared<tf2_ros::TransformBroadcaster>(this);

    // Odometry publisher
    odom_pub_ = create_publisher<nav_msgs::msg::Odometry>("/odom", 10);

    // 20 ms timer
    timer_ = create_wall_timer(20ms, std::bind(&VesselKinematics::update, this));
  }

private:
  void cmdVelCallback(const geometry_msgs::msg::Twist::SharedPtr msg)
  {
    // store commanded body‐frame velocities
    u_ = msg->linear.x;
    v_ = msg->linear.y;
    w_ = msg->linear.z;
    p_ = msg->angular.x;
    q_ = msg->angular.y;
    r_ = msg->angular.z;
  }

  void update()
  {
    const double dt = 0.02;  // 20 ms
    
    // Precompute sines and cosines
    double sr = std::sin(roll_),   cr = std::cos(roll_);
    double sp = std::sin(pitch_),  cp = std::cos(pitch_);
    double sy = std::sin(yaw_),    cy = std::cos(yaw_);
    double tp = std::tan(pitch_);

    // --- 1) Compute r_dot = R_n_b(roll,pitch,yaw) * [u,v,w] ---
    // R_n_b from eqn 12 in the paper
    double x_dot = cy*cp * u_
                 + (cy*sp*sr - sy*cr) * v_
                 + (cy*sp*cr + sy*sr) * w_;
    double y_dot = sy*cp * u_
                 + (sy*sp*sr + cy*cr) * v_
                 + (sy*sp*cr - cy*sr) * w_;
    double z_dot = -sp   * u_
                 + cp*sr * v_
                 + cp*cr * w_;

    // --- 2) Compute Euler derivatives --- 
    // Θ̇ = T_b(Θ) * [p,q,r]   (eqn 26 & discussion)
    // Eq 28: direct Td * [p,q,r]
    double roll_dot  = p_ 
                    + sr*tp * q_
                    + cr*tp * r_;
    double pitch_dot = cr    * q_
                    - sr    * r_;
    double yaw_dot   = sr/cp * q_
                    + cr/cp * r_;

    // --- 3) Integrate to update pose ---
    x_    += x_dot   * dt;
    y_    += y_dot   * dt;
    z_    += z_dot   * dt;
    roll_  = normalize(roll_  + roll_dot  * dt);
    pitch_ = normalize(pitch_ + pitch_dot * dt);
    yaw_   = normalize(yaw_   + yaw_dot   * dt);

    // --- 4) Publish TF map -> base_link ---
    geometry_msgs::msg::TransformStamped tf;
    tf.header.stamp    = now();
    tf.header.frame_id = "map";
    tf.child_frame_id  = "base_link";

    tf.transform.translation.x = x_;
    tf.transform.translation.y = y_;
    tf.transform.translation.z = z_;

    tf2::Quaternion q;
    q.setRPY(roll_, pitch_, yaw_);
    tf.transform.rotation.x = q.x();
    tf.transform.rotation.y = q.y();
    tf.transform.rotation.z = q.z();
    tf.transform.rotation.w = q.w();

    tf_broadcaster_->sendTransform(tf);


    // --- 5) Publish Odometry on /odom ---
    nav_msgs::msg::Odometry odom;
    odom.header.stamp    = tf.header.stamp;
    odom.header.frame_id = "map";
    odom.child_frame_id  = "base_link";

    // pose (inertial)
    odom.pose.pose.position.x    = x_;
    odom.pose.pose.position.y    = y_;
    odom.pose.pose.position.z    = z_;
    odom.pose.pose.orientation.x = q.x();
    odom.pose.pose.orientation.y = q.y();
    odom.pose.pose.orientation.z = q.z();
    odom.pose.pose.orientation.w = q.w();

    // twist (in body frame)
    odom.twist.twist.linear.x  = u_;
    odom.twist.twist.linear.y  = v_;
    odom.twist.twist.linear.z  = w_;
    odom.twist.twist.angular.x = p_;
    odom.twist.twist.angular.y = q_;
    odom.twist.twist.angular.z = r_;

    odom_pub_->publish(odom);
  }

  double normalize(double angle)
  {
    // wrap to [-π, π]
    while (angle > M_PI)  angle -= 2*M_PI;
    while (angle < -M_PI) angle += 2*M_PI;
    return angle;
  }

  // State
  double x_, y_, z_, roll_, pitch_, yaw_;

  // Body‐frame velocity commands
  double u_, v_, w_, p_, q_, r_;

  // ROS interfaces
  rclcpp::Subscription<geometry_msgs::msg::Twist>::SharedPtr sub_cmd_vel_;
  std::shared_ptr<tf2_ros::TransformBroadcaster> tf_broadcaster_;
  rclcpp::Publisher<nav_msgs::msg::Odometry>::SharedPtr odom_pub_;
  rclcpp::TimerBase::SharedPtr timer_;
};

int main(int argc, char **argv)
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<VesselKinematics>());
  rclcpp::shutdown();
  return 0;
}
