// src/vessel_dynamics_node.cpp

#include <memory>
#include <cmath>
#include <rclcpp/rclcpp.hpp>

#include <geometry_msgs/msg/wrench.hpp>
#include <nav_msgs/msg/odometry.hpp>
#include <geometry_msgs/msg/transform_stamped.hpp>
#include <tf2_ros/transform_broadcaster.h>
#include <tf2/LinearMath/Quaternion.h>

#include <Eigen/Dense>

using namespace std::chrono_literals;

/**
 * VesselDynamics integrates:
 *  - Rigid‐body EoM in body frame (paper eq 92):
 *      m·ν̇₁ = m·S(ν₁)·ν₂ + τ_force
 *      I·ν̇₂ =     S(I·ν₂)·ν₂ + τ_moment
 *  - Kinematic map to inertial frame:
 *      ṗ = R(φ,θ,ψ) · ν₁        (paper eq 12)
 *      [φ̇,θ̇,ψ̇]ᵀ = T(φ,θ) · ν₂  (paper eq 28)
 *
 * Outputs both /odom and map→base_link TF.
 */
class VesselDynamics : public rclcpp::Node
{
public:
  VesselDynamics() : Node("vessel_dynamics")
  {
    //----------------------------------------------------------------
    // 1) Parameters
    //----------------------------------------------------------------
    m_ = 500.0; // [kg]
    // principal inertias [kg·m²] from your URDF:
    I_b_.diagonal() << 30.21, 168.54, 181.33;

    //----------------------------------------------------------------
    // 2) State initialization
    //----------------------------------------------------------------
    nu_.setZero(); // [u,v,w,p,q,r]
    // nu_ << 0.0, 0.0, 0.0, 0.0, 0.0, 3.0;
    eta_.setZero(); // [x,y,z]
    // eta_ << 10.0, 0.0, 0.0;
    eul_.setZero(); // [roll, pitch, yaw]
    tau_force_.setZero();
    tau_moment_.setZero();

    //----------------------------------------------------------------
    // 3) ROS interfaces
    //----------------------------------------------------------------
    sub_wrench_ = create_subscription<geometry_msgs::msg::Wrench>(
        "/cmd_force", 10, std::bind(&VesselDynamics::onWrench, this, std::placeholders::_1));

    // up to 10 messages will be buffered if the subscriber can’t keep up
    // we call odom_pub_->publish(odom); in step() to actually send out the odometry data.
    odom_pub_  = create_publisher<nav_msgs::msg::Odometry>("/odom", 10);
    nudot_pub_ = create_publisher<geometry_msgs::msg::Twist>("/nudot", 10); // debug

    // coordinate frame transforms (map → base_link)
    tf_broadcaster_ = std::make_shared<tf2_ros::TransformBroadcaster>(this);

    // calls step() each 20ms
    timer_ = create_wall_timer(20ms, std::bind(&VesselDynamics::step, this));
  }

private:
  // callback: capture incoming body‐frame forces & moments
  void onWrench(const geometry_msgs::msg::Wrench::SharedPtr msg)
  {
    tau_force_ << msg->force.x, msg->force.y, msg->force.z;
    tau_moment_ << msg->torque.x, msg->torque.y, msg->torque.z;
  }

  void step()
  {
    const double dt = 0.02; // 50 Hz

    // Split ν into translational ν₁ and rotational ν₂
    Eigen::Vector3d nu1 = nu_.segment<3>(0); // 3x1 matrix
    Eigen::Vector3d nu2 = nu_.segment<3>(3); // 3x1 matrix

    // ------------------------------------------------------------
    //  A) Rigid‐body dynamics (paper eq 92)
    // ------------------------------------------------------------
    // S(v) = skew-symmetric matrix such that S(v)*w = v×w
    auto S = [&](const Eigen::Vector3d& v)
    {
      Eigen::Matrix3d m;
      m << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
      return m;
    };

    // --- block 1: m·ν̇₁ = m·S(ν₁)·ν₂ + τ_force
    Eigen::Vector3d nu1_dot = (m_ * S(nu1) * nu2 + tau_force_) / m_;

    Eigen::Vector3d I_nu2   = I_b_ * nu2;
    Eigen::Vector3d nu2_dot = I_b_.inverse() * (m_ * S(nu1) * nu1 + S(I_nu2) * nu2 + tau_moment_);

    // assemble ν̇ = [ν̇₁; ν̇₂]
    Eigen::Matrix<double, 6, 1> nudot;
    nudot << nu1_dot, nu2_dot;

    // debug
    geometry_msgs::msg::Twist nudot_msg;
    nudot_msg.linear.x  = nu1_dot.x();
    nudot_msg.linear.y  = nu1_dot.y();
    nudot_msg.linear.z  = nu1_dot.z();
    nudot_msg.angular.x = nu2_dot.x();
    nudot_msg.angular.y = nu2_dot.y();
    nudot_msg.angular.z = nu2_dot.z();
    nudot_pub_->publish(nudot_msg);

    // integrate body‐frame velocity νk+1​=νk​+ν˙k​⋅Δt
    nu_ += nudot * dt;
    nu1 = nu_.segment<3>(0);
    nu2 = nu_.segment<3>(3);

    // ------------------------------------------------------------
    //  B) Kinematic mapping to inertial pose
    // ------------------------------------------------------------
    // 1) Translational
    double cr = std::cos(eul_(0)), sr = std::sin(eul_(0));
    double cp = std::cos(eul_(1)), sp = std::sin(eul_(1));
    double cy = std::cos(eul_(2)), sy = std::sin(eul_(2));

    Eigen::Matrix3d Rnb;
    Rnb << cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr, //
        sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr,    //
        -sp, cp * sr, cp * cr;

    Eigen::Vector3d eta_dot = Rnb * nu1;
    eta_ += eta_dot * dt;

    // 2) Rotational: Θ̇ = T(φ,θ) * ν₂        (paper eq 28)
    Eigen::Matrix3d T;
    T << 1, sr * sp / cp, cr * sp / cp, 0, cr, -sr, 0, sr / cp, cr / cp;

    Eigen::Vector3d eul_dot = T * nu2;
    eul_ += eul_dot * dt;
    normalizeEuler(eul_);

    // ------------------------------------------------------------
    //  C) Publish TF + Odometry
    // ------------------------------------------------------------
    auto stamp = now();

    // 1) TF: map → base_link
    geometry_msgs::msg::TransformStamped tf_msg;
    tf_msg.header.stamp            = stamp;
    tf_msg.header.frame_id         = "map";
    tf_msg.child_frame_id          = "base_link";
    tf_msg.transform.translation.x = eta_(0);
    tf_msg.transform.translation.y = eta_(1);
    tf_msg.transform.translation.z = eta_(2);

    tf2::Quaternion q;
    q.setRPY(eul_(0), eul_(1), eul_(2));
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
    odom.pose.pose.position.x  = eta_(0);
    odom.pose.pose.position.y  = eta_(1);
    odom.pose.pose.position.z  = eta_(2);
    odom.pose.pose.orientation = tf_msg.transform.rotation;

    // twist (body frame)
    odom.twist.twist.linear.x  = nu_(0);
    odom.twist.twist.linear.y  = nu_(1);
    odom.twist.twist.linear.z  = nu_(2);
    odom.twist.twist.angular.x = nu_(3);
    odom.twist.twist.angular.y = nu_(4);
    odom.twist.twist.angular.z = nu_(5);

    odom_pub_->publish(odom);
  }

  // wrap each Euler angle into [-π,π]
  void normalizeEuler(Eigen::Vector3d& e)
  {
    for (int i = 0; i < 3; i++)
    {
      while (e(i) > M_PI)
        e(i) -= 2 * M_PI;
      while (e(i) < -M_PI)
        e(i) += 2 * M_PI;
    }
  }

  // ---------- node state & parameters ------------
  double                      m_;
  Eigen::Matrix3d             I_b_; // principal inertia
  Eigen::Matrix<double, 6, 1> nu_;  // body‐frame twist [u,v,w,p,q,r]
  Eigen::Vector3d             eta_; // inertial pos [x,y,z]
  Eigen::Vector3d             eul_; // Euler angles [φ,θ,ψ]
  Eigen::Vector3d             tau_force_, tau_moment_;

  // ---------- ROS interfaces ---------------------
  rclcpp::Subscription<geometry_msgs::msg::Wrench>::SharedPtr sub_wrench_;
  rclcpp::Publisher<nav_msgs::msg::Odometry>::SharedPtr       odom_pub_;
  rclcpp::Publisher<geometry_msgs::msg::Twist>::SharedPtr     nudot_pub_; // debug
  std::shared_ptr<tf2_ros::TransformBroadcaster>              tf_broadcaster_;
  rclcpp::TimerBase::SharedPtr                                timer_;
};

int main(int argc, char** argv)
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<VesselDynamics>());
  rclcpp::shutdown();
  return 0;
}
