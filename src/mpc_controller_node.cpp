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
#include "vessel_kinematics/mpc_controller.hpp"

using namespace std::chrono_literals; // write time as ms
using namespace casadi;

class MPCNode : public rclcpp::Node
{

public:
  MPCNode() : Node("mpc_controller_node")
  {
    //----------------------------------------------------------------
    // Variable initialization
    //----------------------------------------------------------------
    this->loadParameters();
    mpc_state_.x_k_.setZero();
    mpc_state_.x_k_minus_1_.setZero();
    mpc_state_.tau_k_minus_1_.setZero();
    mpc_state_.A_bar_.setZero();
    mpc_state_.B_bar_.setZero();
    //----------------------------------------------------------------
    // ROS interfaces
    //----------------------------------------------------------------
    pos_cmd_sub_ = create_subscription<geometry_msgs::msg::PoseStamped>(
        "/goal_pose", 10, std::bind(&MPCNode::goalCallback, this, std::placeholders::_1));

    state_sub_ = create_subscription<nav_msgs::msg::Odometry>(
        "/vessel_state", 10, std::bind(&MPCNode::stateCallback, this, std::placeholders::_1));

    wrench_pub_ = create_publisher<geometry_msgs::msg::Wrench>("/cmd_wrench", 10);

    // calls step() each 20ms
    timer_ =
        create_wall_timer(std::chrono::duration<double>(mpc_params_.Ts_), std::bind(&MPCNode::controllerLoop, this));

    // --
    this->setupSolver();
  }

private:
  void goalCallback(const geometry_msgs::msg::PoseStamped::SharedPtr msg)
  {
    // NOTE: X and Y commands of goal_pose comply with ROS convention ENU frame,
    // while yaw is given directly as degrees in NED frame for convenience.
    goal_received_ = true;
    std::cout << goal_received_ << std::endl;
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
    eta_(2) = vk::wrapAngle(M_PI / 2.0 - yaw_enu); // convert ENU yaw (CCW, 0 = East) to NED yaw (CW, 0 = North)
  }
  void controllerLoop()
  {
    if (!goal_received_)
    {
      return;
    }

    // 1. State Preparation
    mpc_state_.x_k_ << nu_, eta_;
    Eigen::Vector<double, 6> delta_x_k_ = mpc_state_.x_k_ - mpc_state_.x_k_minus_1_;
    Eigen::Vector<double, 9> p_k_;
    p_k_ << delta_x_k_, eta_;

    // 2. Model Calculation
    // 2.1 continuous matrices
    Eigen::Matrix<double, 6, 6> A_;
    A_.setZero();
    vk::updateCoriolisMatrix(nu_, this->C_);
    vk::updateDampingMatrix(nu_, this->D_);
    A_.topLeftCorner<3, 3>()     = -vk::M_.ldlt().solve(C_ + D_); // efficient way to calculate
    A_.bottomLeftCorner<3, 3>()  = vk::calculateRotationMatrix(eta_);
    A_.bottomRightCorner<3, 3>() = vk::calculateKinematicJacobian(nu_, eta_);

    Eigen::Matrix<double, 6, 3> B_;
    B_.setZero();
    B_.topLeftCorner<3, 3>() = vk::M_.inverse();

    // 2.2 discretized matrices
    Eigen::Matrix<double, 6, 6> A_d_ = Eigen::Matrix<double, 6, 6>::Identity() + mpc_params_.Ts_ * A_;
    Eigen::Matrix<double, 6, 3> B_d_ = mpc_params_.Ts_ * B_;
    Eigen::Matrix<double, 3, 6> C_d_;
    C_d_ << Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Identity();

    // 2.3 augmented matrices
    mpc_state_.A_bar_.topLeftCorner<6, 6>()     = A_d_;
    mpc_state_.A_bar_.bottomLeftCorner<3, 6>()  = C_d_ * A_d_;
    mpc_state_.A_bar_.bottomRightCorner<3, 3>() = Eigen::Matrix3d::Identity();

    mpc_state_.B_bar_.topRows<6>()    = B_d_;
    mpc_state_.B_bar_.bottomRows<3>() = C_d_ * B_d_;

    // 3. Set Solver Parameters
    // 3.1 Conversion
    std::vector<double> p_k_vec(p_k_.data(),
                                p_k_.data() + p_k_.size()); // TODO learn the typed pointed scaled addition logic
    casadi::DM          p_k_dm = casadi::DM(p_k_vec);

    Eigen::Vector<double, 9> p_ref;
    p_ref.setZero();
    p_ref.tail<3>() = eta_des_;
    std::vector<double> p_ref_vec(p_ref.data(), p_ref.data() + p_ref.size());
    casadi::DM          p_ref_dm = casadi::DM(p_ref_vec);

    std::vector<double> tau_prev_vec(mpc_state_.tau_k_minus_1_.data(),
                                     mpc_state_.tau_k_minus_1_.data() + mpc_state_.tau_k_minus_1_.size());
    casadi::DM          tau_prev_dm = casadi::DM(tau_prev_vec);

    std::vector<double> A_bar_vec(mpc_state_.A_bar_.data(), mpc_state_.A_bar_.data() + mpc_state_.A_bar_.size());
    casadi::DM          A_bar_dm = casadi::DM::reshape(casadi::DM(A_bar_vec), 9, 9);

    std::vector<double> B_bar_vec(mpc_state_.B_bar_.data(), mpc_state_.B_bar_.data() + mpc_state_.B_bar_.size());
    casadi::DM          B_bar_dm = casadi::DM::reshape(casadi::DM(B_bar_vec), 9, 3);

    // 3.2 Send to the optimizer
    mpc_solver_.opti_->set_value(mpc_solver_.p0_param_, p_k_dm);
    mpc_solver_.opti_->set_value(mpc_solver_.p_ref_param_, p_ref_dm);
    mpc_solver_.opti_->set_value(mpc_solver_.tau_prev_param_, tau_prev_dm);
    mpc_solver_.opti_->set_value(mpc_solver_.A_bar_param_, A_bar_dm);
    mpc_solver_.opti_->set_value(mpc_solver_.B_bar_param_, B_bar_dm);

    // 4. Send to the solver
    casadi::OptiSol solution = mpc_solver_.opti_->solve();

    // 5. Extract and publish
    casadi::DM          delta_tau_k_dm  = solution.value(mpc_solver_.U_var_)(Slice(), 0);
    std::vector<double> delta_tau_k_vec = (std::vector<double>)delta_tau_k_dm;
    Eigen::Vector3d     delta_tau_k     = Eigen::Map<Eigen::Vector3d>(delta_tau_k_vec.data());

    // Calculate the new absolute command
    Eigen::Vector3d tau_k = mpc_state_.tau_k_minus_1_ + delta_tau_k;

    geometry_msgs::msg::Wrench cmd;
    cmd.force.x  = tau_k(0); // equivalent to tau_des_(0)
    cmd.force.y  = tau_k(1); // equivalent to tau_des_(0)
    cmd.torque.z = tau_k(2); // equivalent to tau_des_(0)
    wrench_pub_->publish(cmd);

    // 6. Update state for next loop
    mpc_state_.x_k_minus_1_   = mpc_state_.x_k_;
    mpc_state_.tau_k_minus_1_ = tau_k;
  }
  void setupSolver()
  {
    mpc_solver_.opti_ = std::make_shared<casadi::Opti>();
    const int p_dim   = 9; // Dimension of augmented state p = [Δx, y]
    const int u_dim   = 3; // Dimension of control input u = Δτ

    // Variables the solver will find
    mpc_solver_.U_var_ = mpc_solver_.opti_->variable(u_dim, mpc_params_.P_);
    mpc_solver_.P_var_ = mpc_solver_.opti_->variable(p_dim, mpc_params_.P_ + 1);

    // Parameters we will provide at runtime
    mpc_solver_.p0_param_       = mpc_solver_.opti_->parameter(p_dim, 1);
    mpc_solver_.p_ref_param_    = mpc_solver_.opti_->parameter(p_dim, 1);
    mpc_solver_.tau_prev_param_ = mpc_solver_.opti_->parameter(u_dim, 1);
    mpc_solver_.A_bar_param_    = mpc_solver_.opti_->parameter(p_dim, p_dim);
    mpc_solver_.B_bar_param_    = mpc_solver_.opti_->parameter(p_dim, u_dim);

    // --- Symbolic Cost Function ---
    casadi::MX cost = 0;
    for (int i = 0; i < mpc_params_.P_; i++)
    {
      auto state_error    = mpc_solver_.P_var_(Slice(), i + 1) - mpc_solver_.p_ref_param_;
      auto control_change = mpc_solver_.U_var_(Slice(), i);
      cost += casadi::MX::mtimes({state_error.T(), mpc_params_.Q_, state_error});
      cost += casadi::MX::mtimes({control_change.T(), mpc_params_.R_, control_change});
    }
    mpc_solver_.opti_->minimize(cost);

    // --- Define Prediction Constraints ---
    // Initial state constraint
    mpc_solver_.opti_->subject_to(mpc_solver_.P_var_(Slice(), 0) == mpc_solver_.p0_param_);
    // Dynamic constraints over the horizon
    // p[k+i+1] = Ā * p[k+i] + B̄ * Δτ[k+i]
    for (int i = 0; i < mpc_params_.P_; i++)
    {
      mpc_solver_.opti_->subject_to(mpc_solver_.P_var_(Slice(), i + 1) ==
                                    casadi::MX::mtimes(mpc_solver_.A_bar_param_, mpc_solver_.P_var_(Slice(), i)) +
                                        casadi::MX::mtimes(mpc_solver_.B_bar_param_, mpc_solver_.U_var_(Slice(), i)));
    }
    // --- Define Control Constraints ---
    casadi::MX tau_i = mpc_solver_.tau_prev_param_; // last known thrust
    for (int i = 0; i < mpc_params_.P_; i++)
    {
      tau_i += mpc_solver_.U_var_(Slice(), i);

      mpc_solver_.opti_->subject_to(tau_i >= mpc_params_.tau_min_);
      mpc_solver_.opti_->subject_to(tau_i <= mpc_params_.tau_max_);
    }
    // --- Create the Solver ---
    casadi::Dict opts;
    opts["qpsol"]         = "qpoases";
    opts["qpsol_options"] = casadi::Dict{{"printLevel", "none"}};
    opts["print_time"]    = false;
    opts["expand"]        = true;

    mpc_solver_.opti_->solver("sqpmethod", opts);
    RCLCPP_INFO(this->get_logger(), "MPC solver has been set up successfully.");
  }
  void loadParameters()
  {
    MPCParams p_default;

    mpc_params_.P_  = this->declare_parameter<int>("P", p_default.P_);
    mpc_params_.Ts_ = this->declare_parameter<double>("Ts", p_default.Ts_);

    std::vector<double> Q_vec       = this->get_parameter("Q_weights").as_double_array();
    std::vector<double> R_vec       = this->get_parameter("R_weights").as_double_array();
    std::vector<double> tau_min_vec = this->get_parameter("tau_min").as_double_array();
    std::vector<double> tau_max_vec = this->get_parameter("tau_max").as_double_array();
    mpc_params_.Q_                  = casadi::DM::diag(Q_vec);
    mpc_params_.R_                  = casadi::DM::diag(R_vec);
    mpc_params_.tau_min_            = casadi::DM::reshape(tau_min_vec, 3, 1);
    mpc_params_.tau_max_            = casadi::DM::reshape(tau_max_vec, 3, 1);
  }
  // ROS2 Subscribers and Publishers
  rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr pos_cmd_sub_;
  rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr         state_sub_;
  rclcpp::Publisher<geometry_msgs::msg::Wrench>::SharedPtr         wrench_pub_;
  rclcpp::TimerBase::SharedPtr                                     timer_;

  // Variable declarations
  Eigen::Vector3d eta_;     // position [x, y, ψ]
  Eigen::Vector3d eta_des_; // [N, E, ψ]
  Eigen::Vector3d nu_;      // body-frame twist [u, v, r]
  Eigen::Matrix3d C_;
  Eigen::Matrix3d D_;
  MPCParams       mpc_params_;
  MPCState        mpc_state_;
  MPCSolver       mpc_solver_;

  bool goal_received_ = false;
};

int main(int argc, char** argv)
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<MPCNode>());
  rclcpp::shutdown();
  return 0;
}
