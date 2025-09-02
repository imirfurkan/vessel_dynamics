#include <memory>
#include <chrono>
#include <rclcpp/rclcpp.hpp>
#include <nav_msgs/msg/odometry.hpp>
#include "geometry_msgs/msg/wrench.hpp"
#include <geometry_msgs/msg/pose_stamped.hpp>
#include "std_msgs/msg/float64_multi_array.hpp"
#include "vessel_kinematics/thrust_allocator.hpp"
#include "vessel_kinematics/utils.hpp"
#include <tf2/LinearMath/Quaternion.h>
#include <tf2/LinearMath/Matrix3x3.h>

class ControllerNode : public rclcpp::Node
{

public:
  ControllerNode() : Node("milliampere_controller_node")
  {
    //----------------------------------------------------------------
    // State initialization
    //----------------------------------------------------------------
    nu_.setZero();
    eta_dot_ref_.setZero();
    eta_ddot_ref_.setZero();

    Kp.setZero();
    Ki.setZero();
    Kd.setZero();
    Kp.diagonal() << 200, 200, 800;
    Ki.diagonal() << 15, 15, 45;
    Kd.diagonal() << 700, 700, 1600;

    // Kp.diagonal() << 0, 0, 170;
    // Ki.diagonal() << 0, 0, 4;
    // Kd.diagonal() << 0, 0, 1600;

    // Kp.diagonal() << 200, 200, 260;
    // Ki.diagonal() << 0, 0, 0;
    // Kd.diagonal() << 700, 700, 1800;

    integral_error_limits_ << 10.0, 10.0, 1.0; // TODO tune these

    // TODO check if the logic & limits makes sense, in the way that if reference model function always gives the max
    // values written here.
    vel_ref_max_ << 2.57,    // Max speed ~ 5 knots [m/s]
        1.0,                 // TODO (arbitrary)
        32.0 * M_PI / 180.0; // Table 5.1 (not azimuth but body?)

    acc_ref_max_ << 0.25, // TODO (arbitrary)
        0.1,              // TODO (arbitrary)
        0.2;              // TODO (arbitrary)

    // tuning parameters for the reference model
    // reference model should be a slow, smooth system that ignores fast disturbances
    // as we want reference feedforward model to be low bandwidth and pid to be high bandwidth
    ref_zeta_ << 1.0, 1.0, 1.0;     // Critically damped is usually best
    ref_omega_n_ << 0.5, 0.25, 0.5; // Tune these for response speed, higher values create
                                    // high-bandwidth reference signal

    // // --- Time ---
    // last_time_ = this->get_clock()->now();

    //----------------------------------------------------------------
    // ROS interfaces
    //----------------------------------------------------------------
    pos_cmd_sub_ = create_subscription<geometry_msgs::msg::PoseStamped>(
        "/goal_pose", 10, std::bind(&ControllerNode::goalCallback, this, std::placeholders::_1)); // TODO anla

    state_sub_ = create_subscription<nav_msgs::msg::Odometry>(
        "/vessel_state", 10, std::bind(&ControllerNode::stateCallback, this, std::placeholders::_1));

    wrench_pub_ = create_publisher<geometry_msgs::msg::Wrench>("/cmd_wrench", 10);
  }

private:
  // callback: capture incoming goal position
  // TODO learn how multiple callbacks are handled, in which order and frequency etc.
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

    // TODO without implementing a logic like the below, we shouldn't use the program for multiple different goal
    // positions without reset
    // // A new goal has just been received. Reset the reference model to the current state.
    // RCLCPP_INFO(this->get_logger(), "New goal received! Resetting reference model to current vessel state.");
    // eta_ref_ = eta_;
    // eta_dot_ref_.setZero();
    // eta_ddot_ref_.setZero();
  }

  void stateCallback(const nav_msgs::msg::Odometry::SharedPtr msg)
  {
    // Dynamic Time calculation. After test, I saw that it's ok to just use vk::dt as dynamic time is also very similar.
    // rclcpp::Time current_time = this->get_clock()->now();
    // double       dt           = (current_time - last_time_).seconds();
    // last_time_                = current_time;
    // RCLCPP_INFO(this->get_logger(), "dynamic time: [%f]", dt);

    // // On the first run, dt can be large and meaningless, so we skip the control logic
    // if (first_run_)
    // {
    //   first_run_ = false;
    //   return;
    // }

    // Odometry Twist
    nu_(0) = msg->twist.twist.linear.x;  // u
    nu_(1) = msg->twist.twist.linear.y;  // v
    nu_(2) = msg->twist.twist.angular.z; // r
    // Odometry Pose is published in ENU frame, convert to NED
    eta_(0) = msg->pose.pose.position.y;
    eta_(1) = msg->pose.pose.position.x;

    tf2::Quaternion q(msg->pose.pose.orientation.x,
                      msg->pose.pose.orientation.y,
                      msg->pose.pose.orientation.z,
                      msg->pose.pose.orientation.w);

    double roll, pitch, yaw_enu;
    tf2::Matrix3x3(q).getRPY(roll, pitch, yaw_enu);

    // convert ENU yaw (CCW, 0 = East) to NED yaw (CW, 0 = North)
    eta_(2) = vk::wrapAngle(M_PI / 2.0 - yaw_enu);

    if (goal_received_)
    {
      // Record the start time using std::chrono::steady_clock
      auto start_time = std::chrono::steady_clock::now();
      // ---------- CONTROL PIPELINE
      Eigen::Matrix3d Rbn = vk::calculateRotationMatrix(eta_).transpose(); // of eta_, not eta_des_

      updateReferenceModel(vk::dt);
      nu_ref_     = Rbn * eta_dot_ref_;  // velocity in the body frame
      nu_dot_ref_ = Rbn * eta_ddot_ref_; // acceleration in the body frame
      // RCLCPP_INFO(this->get_logger(), "nu_ref_: [%f, %f, %f]", nu_ref_(0), nu_ref_(1), nu_ref_(2));
      // RCLCPP_INFO(this->get_logger(), "nu_dot_ref_: [%f, %f, %f]", nu_dot_ref_(0), nu_dot_ref_(1), nu_dot_ref_(2));

      Eigen::Vector3d tau_ff  = calculateFF();
      Eigen::Vector3d tau_pid = calculatePID(Rbn, vk::dt);

      tau_des_ = tau_pid + tau_ff;

      // RCLCPP_INFO(this->get_logger(), "FF Tau: [%f, %f, %f]", tau_ff(0), tau_ff(1), tau_ff(2));
      // RCLCPP_INFO(this->get_logger(), "PID Tau: [%f, %f, %f]", tau_pid(0), tau_pid(1), tau_pid(2));
      // RCLCPP_INFO(this->get_logger(), "Desired Tau: [%f, %f, %f]", tau_des_(0), tau_des_(1), tau_des_(2));

      geometry_msgs::msg::Wrench cmd;
      cmd.force.x  = tau_des_(0);
      cmd.force.y  = tau_des_(1);
      cmd.torque.z = tau_des_(2);
      wrench_pub_->publish(cmd);

      // Record the end time
      auto end_time = std::chrono::steady_clock::now();

      // Calculate the duration
      std::chrono::duration<double, std::milli> duration = end_time - start_time;

      // Print the elapsed time
      std::cout << "Controller took: " << duration.count() << " milliseconds." << std::endl;
    }
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
    eta_ddot_ref_ = vk::clampVec3(eta_ddot_ref_, -acc_ref_max_, acc_ref_max_); // Saturate the acceleration
    eta_dot_ref_ += eta_ddot_ref_ * dt;
    eta_dot_ref_ = vk::clampVec3(eta_dot_ref_, -vel_ref_max_, vel_ref_max_);
    eta_ref_ += eta_dot_ref_ * dt;
    eta_ref_(2) = vk::wrapAngle(eta_ref_(2));
  }

  Eigen::Vector3d calculatePID(Eigen::Matrix3d Rbn, double dt)
  {
    Eigen::Vector3d eta_diff_;
    eta_diff_ << (eta_(0) - eta_ref_(0)),     // N error
        (eta_(1) - eta_ref_(1)),              // E error
        vk::wrapAngle(eta_(2) - eta_ref_(2)); // wrapped ψ error

    Eigen::Vector3d error_     = Rbn * eta_diff_; // error should be in the body frame
    Eigen::Vector3d error_dot_ = nu_ - nu_ref_;

    error_integral_ += error_ * dt;
    std::cout << "Error integral before clamping: " << error_integral_(2) << std::endl;
    error_integral_ = vk::clampVec3(error_integral_, -integral_error_limits_, integral_error_limits_);

    Eigen::Vector3d tau = -Kp * error_ - Ki * error_integral_ - Kd * error_dot_;
    return tau;
  }

  Eigen::Vector3d calculateFF()
  {
    Eigen::Matrix3d C_d, D_d;
    vk::updateCoriolisMatrix(nu_ref_, C_d);
    vk::updateDampingMatrix(nu_ref_, D_d);
    Eigen::Vector3d tau = vk::M_ * nu_dot_ref_ + C_d * nu_ref_ + D_d * nu_ref_;
    return tau;
  }
  // --- Time ---
  rclcpp::Time last_time_;
  bool         first_run_ = true;

  // ROS2 Subscribers and Publishers
  rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr pos_cmd_sub_;
  rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr         state_sub_;
  rclcpp::Publisher<geometry_msgs::msg::Wrench>::SharedPtr         wrench_pub_;

  // Variable declarations // TODO check variable declarations, needed here or there etc
  // TODO when to make them member variables vs
  Eigen::Vector3d eta_;                   // position [x, y, ψ]
  Eigen::Vector3d vel_ref_max_;           // v_max in the book // physical limit
  Eigen::Vector3d acc_ref_max_;           // a_max (derived from ν̇_max) // physical limit
  Eigen::Vector3d ref_omega_n_;           // Diagonal elements of Ω
  Eigen::Vector3d ref_zeta_;              // Diagonal elements of Δ
  Eigen::Vector3d eta_ref_;               // η_d   (smooth desired position)
  Eigen::Vector3d eta_dot_ref_;           // η̇_d or v_d (smooth desired velocity)
  Eigen::Vector3d eta_ddot_ref_;          // η̈_d or a_d (smooth desired acceleration)
  Eigen::Vector3d eta_des_;               // [N, E, ψ]
  Eigen::Vector3d tau_des_;               // [X, Y, N]
  Eigen::Matrix3d Kp;                     //
  Eigen::Matrix3d Ki;                     //
  Eigen::Matrix3d Kd;                     //
  Eigen::Vector3d integral_error_limits_; //
  Eigen::Vector3d error_integral_{};      //
  Eigen::Vector3d nu_;                    // body-frame twist [u, v, r]
  Eigen::Vector3d nu_ref_;
  Eigen::Vector3d nu_dot_ref_;

  bool goal_received_ = false;
};

int main(int argc, char** argv)
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<ControllerNode>());
  rclcpp::shutdown();
  return 0;
}
