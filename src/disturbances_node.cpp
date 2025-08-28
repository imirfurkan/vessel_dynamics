#include <rclcpp/rclcpp.hpp>
#include <geometry_msgs/msg/wrench.hpp>
#include <Eigen/Dense>
#include <random>
#include <chrono>
#include "vessel_kinematics/utils.hpp"

using namespace std::chrono_literals; // write time as ms

class DisturbancesNode : public rclcpp::Node
{
public:
  DisturbancesNode() : Node("disturbances_node"), distribution_(0.0, 1.0) // Initialize distribution here
  // TODO why is distribution_ here
  {
    //----------------------------------------------------------------
    // Wave State initialization
    //----------------------------------------------------------------
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
    // ROS interfaces
    //----------------------------------------------------------------
    disturbance_pub_ = this->create_publisher<geometry_msgs::msg::Wrench>("/wave_wrench", 10);
    timer_           = this->create_wall_timer(20ms, std::bind(&DisturbancesNode::timerCallback, this));
    // Seed the random number generator ONCE in the constructor.
    generator_.seed(std::chrono::system_clock::now().time_since_epoch().count());
  }

private:
  void timerCallback()
  {
    std::default_random_engine       generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::normal_distribution<double> distribution(0.0, 1.0); // Mean 0, StdDev 1

    Eigen::Vector3d white_noise_ = {distribution(generator), distribution(generator), distribution(generator)};

    Eigen::Vector<double, 6> xi_dot = A_wave_ * xi_wave_ + E_wave_ * white_noise_;
    xi_wave_ += xi_dot * vk::dt;

    Eigen::Vector3d first_order_waves = K_wave_.cwiseProduct(xi_wave_.tail<3>());
    tau_wave_                         = first_order_waves; // TODO + second_order_waves

    auto wrench_msg     = geometry_msgs::msg::Wrench();
    wrench_msg.force.x  = tau_wave_(0);
    wrench_msg.force.y  = tau_wave_(1);
    wrench_msg.force.z  = 0.0;
    wrench_msg.torque.x = 0.0;
    wrench_msg.torque.y = 0.0;
    wrench_msg.torque.z = tau_wave_(2);
    disturbance_pub_->publish(wrench_msg); // TODO understand the pointer
  }

  // --- ROS 2 Member Variables ---
  rclcpp::Publisher<geometry_msgs::msg::Wrench>::SharedPtr disturbance_pub_;
  rclcpp::TimerBase::SharedPtr                             timer_;

  // ---------- Wave Disturbance Model ------------
  std::default_random_engine       generator_;
  std::normal_distribution<double> distribution_;
  Eigen::Vector3d                  omega0_;   // Dominant wave frequencies [rad/s] for surge, sway, yaw
  Eigen::Vector3d                  lambda_;   // Damping ratios for surge, sway, yaw
  Eigen::Vector3d                  K_wave_;   // Gains to control force/torque magnitude
  Eigen::Vector<double, 6>         xi_wave_;  // Grouped state vector [pos_s, pos_y, pos_n, vel_s, vel_y, vel_n]
  Eigen::Matrix<double, 6, 6>      A_wave_;   // 6x6 state matrix
  Eigen::Matrix<double, 6, 3>      E_wave_;   // 6x3 input matrix
  Eigen::Vector3d                  tau_wave_; // [X, Y, N]
};

int main(int argc, char** argv)
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<DisturbancesNode>());
  rclcpp::shutdown();
  return 0;
}
