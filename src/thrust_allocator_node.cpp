#include <memory>
#include <rclcpp/rclcpp.hpp>
#include "geometry_msgs/msg/wrench.hpp"
#include "std_msgs/msg/float64_multi_array.hpp"
#include "vessel_kinematics/thrust_allocator.hpp"
#include "vessel_kinematics/utils.hpp"

class ThrustAllocatorNode : public rclcpp::Node
{

public:
  ThrustAllocatorNode() : Node("thrust_allocator_node")
  {
    this->loadParameters(); // TODO understand this->

    wrench_sub_ = this->create_subscription<geometry_msgs::msg::Wrench>(
        "/cmd_wrench",
        10,
        std::bind(&ThrustAllocatorNode::wrenchCallBack, this, std::placeholders::_1)); // TODO understand

    actual_wrench_pub_ = this->create_publisher<geometry_msgs::msg::Wrench>("/actual_wrench", 10);
    thruster_cmds_pub_ = this->create_publisher<std_msgs::msg::Float64MultiArray>("/thruster_commands", 10);
  }

private:
  void loadParameters()
  { // Use a temporary default params struct to get default values for declarations
    TAParams p_default;

    ta_params_.Lx = this->declare_parameter<double>("Lx", p_default.Lx);

    // For vector parameters, we get them as std::vector<double> and convert to Eigen
    // std::vector<double> f_min_vec, f_max_vec, dalpha_min_vec, dalpha_max_vec, Wf_vec, Qs_vec, Oa_vec;

    auto f_min_vec = this->declare_parameter<std::vector<double>>(
        "f_min", {p_default.f_min(0), p_default.f_min(1)});           // TODO understand
    ta_params_.f_min = Eigen::Map<Eigen::Vector2d>(f_min_vec.data()); // TODO understand .data(), pointers, copying etc.

    auto f_max_vec   = this->declare_parameter<std::vector<double>>("f_max", {p_default.f_max(0), p_default.f_max(1)});
    ta_params_.f_max = Eigen::Map<Eigen::Vector2d>(f_max_vec.data());

    auto dalpha_min_vec = this->declare_parameter<std::vector<double>>(
        "dalpha_min_rate_rad_s", {p_default.dalpha_min_rate_rad_s(0), p_default.dalpha_min_rate_rad_s(1)});
    ta_params_.dalpha_min_rate_rad_s = Eigen::Map<Eigen::Vector2d>(dalpha_min_vec.data());

    auto dalpha_max_vec = this->declare_parameter<std::vector<double>>(
        "dalpha_max_rate_rad_s", {p_default.dalpha_max_rate_rad_s(0), p_default.dalpha_max_rate_rad_s(1)});
    ta_params_.dalpha_max_rate_rad_s = Eigen::Map<Eigen::Vector2d>(dalpha_max_vec.data());

    auto Wf_vec   = this->declare_parameter<std::vector<double>>("Wf", {p_default.Wf(0), p_default.Wf(1)});
    ta_params_.Wf = Eigen::Map<Eigen::Vector2d>(Wf_vec.data());

    auto Qs_vec =
        this->declare_parameter<std::vector<double>>("Qs", {p_default.Qs(0), p_default.Qs(1), p_default.Qs(2)});
    ta_params_.Qs = Eigen::Map<Eigen::Vector3d>(Qs_vec.data());

    auto Oa_vec   = this->declare_parameter<std::vector<double>>("Oa", {p_default.Oa(0), p_default.Oa(1)});
    ta_params_.Oa = Eigen::Map<Eigen::Vector2d>(Oa_vec.data());
  }

  void wrenchCallBack(const geometry_msgs::msg::Wrench::SharedPtr msg) // TODO understand smart pointer SharedPtr
  {
    Eigen::Vector3d tau_desired;
    tau_desired(0) = msg->force.x;
    tau_desired(1) = msg->force.y;
    tau_desired(2) = msg->torque.z;

    TAResult result = allocate_tau(tau_desired, ta_state_, ta_params_, vk::dt);

    if (result.success)
    {
      ta_state_.alpha = result.alpha;
      ta_state_.f     = result.f;
    }
    else
    {
      ta_state_.f.setZero();
      ta_state_.alpha.setZero();
      RCLCPP_WARN(this->get_logger(), "Thrust allocation failed! Commanding zero thrust and azimuth angles.");
    }

    publishResults();
  }

  void publishResults()
  {
    auto thruster_cmds_msg = std_msgs::msg::Float64MultiArray();
    thruster_cmds_msg.data.resize(4); // TODO understand .data and the need for resize, what was it initially
    thruster_cmds_msg.data[0] = ta_state_.alpha(0); // Azimuth 1 [rad]
    thruster_cmds_msg.data[1] = ta_state_.alpha(1); // Azimuth 2 [rad]
    thruster_cmds_msg.data[2] = ta_state_.f(0);     // Thrust 1 [N]
    thruster_cmds_msg.data[3] = ta_state_.f(1);     // Thrust 2 [N]
    thruster_cmds_pub_->publish(thruster_cmds_msg); // TODO understand the pointer

    auto                        actual_wrench_msg = geometry_msgs::msg::Wrench();
    Eigen::Matrix<double, 3, 2> T                 = Tmatrix(ta_state_.alpha, ta_params_.Lx);
    Eigen::Vector3d             tau_actual        = T * ta_state_.f;
    actual_wrench_msg.force.x                     = tau_actual(0);
    actual_wrench_msg.force.y                     = tau_actual(1);
    actual_wrench_msg.torque.z                    = tau_actual(2);
    actual_wrench_pub_->publish(actual_wrench_msg);

    // debug
    RCLCPP_INFO(this->get_logger(),
                "Azimuths: [%.2f, %.2f], Thrusts: [%.2f, %.2f]",
                ta_state_.alpha(0) * 180 / M_PI,
                ta_state_.alpha(1) * 180 / M_PI,
                ta_state_.f(0),
                ta_state_.f(1));

    RCLCPP_INFO(
        this->get_logger(), "Actual Tau: [X=%.2f, Y=%.2f, N=%.2f]", tau_actual(0), tau_actual(1), tau_actual(2));
  }

  // ROS2 Subscribers and Publishers
  rclcpp::Subscription<geometry_msgs::msg::Wrench>::SharedPtr    wrench_sub_;
  rclcpp::Publisher<geometry_msgs::msg::Wrench>::SharedPtr       actual_wrench_pub_;
  rclcpp::Publisher<std_msgs::msg::Float64MultiArray>::SharedPtr thruster_cmds_pub_;

  // Variables for the allocation algorithm
  TAState  ta_state_;
  TAParams ta_params_; // TODO what is different between this and p_default if both are initialized from the header
                       // file, why did we have to do allat?
};

int main(int argc, char** argv) // TODO understand function parameters and the whole int main function
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<ThrustAllocatorNode>());
  rclcpp::shutdown();
  return 0;
}
