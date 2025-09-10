#include <memory>
#include <rclcpp/rclcpp.hpp>
#include <nav_msgs/msg/odometry.hpp>
#include "geometry_msgs/msg/wrench.hpp"
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <nav_msgs/msg/path.hpp>
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
    RCLCPP_INFO(this->get_logger(), "Starting constructor...");
    this->loadParameters();
    RCLCPP_INFO(this->get_logger(), "Parameters loaded. Initializing state...");
    mpc_state_.x_k_.setZero();
    mpc_state_.x_k_minus_1_.setZero();
    mpc_state_.tau_k_minus_1_.setZero();
    mpc_state_.A_bar_.setZero();
    mpc_state_.B_bar_.setZero();
    RCLCPP_INFO(this->get_logger(), "Parameters loaded. Initializing state...");
    //----------------------------------------------------------------
    // ROS interfaces
    //----------------------------------------------------------------
    pos_cmd_sub_ = create_subscription<geometry_msgs::msg::PoseStamped>(
        "/goal_pose", 10, std::bind(&MPCNode::goalCallback, this, std::placeholders::_1));

    state_sub_ = create_subscription<nav_msgs::msg::Odometry>(
        "/vessel_state", 10, std::bind(&MPCNode::stateCallback, this, std::placeholders::_1));

    wrench_pub_ = create_publisher<geometry_msgs::msg::Wrench>("/cmd_wrench", 10); // TODO revert to cmd_wrench

    predicted_path_pub_ = this->create_publisher<nav_msgs::msg::Path>("/mpc_predicted_path", 10);
    prediction_pub_     = this->create_publisher<nav_msgs::msg::Odometry>("/mpc_prediction", 10);
    // calls step() each 20ms
    timer_ =
        create_wall_timer(std::chrono::duration<double>(mpc_params_.Ts_), std::bind(&MPCNode::controllerLoop, this));
    RCLCPP_INFO(this->get_logger(), "ROS interfaces created. Setting up solver...");
    // --
    this->setupSolver();
    RCLCPP_INFO(this->get_logger(), "Constructor finished successfully!");
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

    // After calculating continuous system matrices
    static bool systems_analyzed = false;
    if (!systems_analyzed)
    {
      // 1. Analyze continuous system (A, B)
      std::string continuous_analysis = vk::getSystemAnalysis(A_, B_, C_d_);
      RCLCPP_INFO(this->get_logger(), "\nContinuous System Analysis:\n%s", continuous_analysis.c_str());

      // 2. Analyze discrete system (A_d, B_d)
      std::string discrete_analysis = vk::getSystemAnalysis(A_d_, B_d_, C_d_);
      RCLCPP_INFO(this->get_logger(), "\nDiscrete System Analysis:\n%s", discrete_analysis.c_str());

      // 3. Analyze augmented system (A_bar, B_bar)
      mpc_state_.A_bar_.topLeftCorner<6, 6>()     = A_d_;
      mpc_state_.A_bar_.bottomLeftCorner<3, 6>()  = C_d_ * A_d_;
      mpc_state_.A_bar_.bottomRightCorner<3, 3>() = Eigen::Matrix3d::Identity();

      mpc_state_.B_bar_.topRows<6>()    = B_d_;
      mpc_state_.B_bar_.bottomRows<3>() = C_d_ * B_d_;

      // Create identity matrix for augmented output matrix
      Eigen::MatrixXd C_bar = Eigen::MatrixXd::Identity(9, 9);

      std::string augmented_analysis = vk::getSystemAnalysis(mpc_state_.A_bar_, mpc_state_.B_bar_, C_bar);
      RCLCPP_INFO(this->get_logger(), "\nAugmented System Analysis:\n%s", augmented_analysis.c_str());

      // Print warnings for critical properties
      if (!vk::isControllable(mpc_state_.A_bar_, mpc_state_.B_bar_))
      {
        RCLCPP_ERROR(this->get_logger(), "Augmented system is not controllable!");
      }
      if (!vk::isObservable(mpc_state_.A_bar_, C_bar))
      {
        RCLCPP_ERROR(this->get_logger(), "Augmented system is not observable!");
      }
      if (!vk::isStable(mpc_state_.A_bar_))
      {
        RCLCPP_WARN(this->get_logger(), "Augmented system contains integrators - marginal stability is expected");
      }

      systems_analyzed = true;
    }

    // Optional: Add periodic stability monitoring
    static rclcpp::Time last_check = this->now();
    if ((this->now() - last_check).seconds() > 5.0)
    { // Check every 5 seconds
      if (!vk::isStable(mpc_state_.A_bar_))
      {
        RCLCPP_WARN_THROTTLE(this->get_logger(),
                             *this->get_clock(),
                             5000, // Throttle to once per 5 seconds
                             "Runtime check: Augmented system stability requires attention");
      }
      last_check = this->now();
    }

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

    //////////////////////////////////////////////////////////////
    // --- NEW: Publish the predicted trajectory for debugging ---
    //////////////////////////////////////////////////////////////

    // Extract the full predicted state sequence from the solver
    casadi::DM predicted_p_sequence = solution.value(mpc_solver_.P_var_);

    // Create a Path message
    auto path_msg            = nav_msgs::msg::Path();
    path_msg.header.stamp    = this->get_clock()->now();
    path_msg.header.frame_id = "map"; // Assuming your vessel state is in the map frame

    // Loop through the prediction horizon
    for (int i = 0; i < mpc_params_.P_ + 1; ++i)
    {
      // Extract the state at step 'i' from the CasADi matrix
      casadi::DM          p_i     = predicted_p_sequence(Slice(), i);
      std::vector<double> p_i_vec = (std::vector<double>)p_i;

      // The position eta is the last 3 elements of the 9D vector p
      double ned_x   = p_i_vec[6];
      double ned_y   = p_i_vec[7];
      double ned_psi = p_i_vec[8];

      // Create a PoseStamped message for this point in the path
      auto pose   = geometry_msgs::msg::PoseStamped();
      pose.header = path_msg.header;

      // Convert from NED (your model's frame) to ENU (RViz's frame)
      pose.pose.position.x = ned_y; // ENU x = East = NED y
      pose.pose.position.y = ned_x; // ENU y = North = NED x

      tf2::Quaternion q;
      double          yaw_enu = vk::wrapAngle(M_PI / 2.0 - ned_psi);
      q.setRPY(0, 0, yaw_enu);
      pose.pose.orientation.x = q.x();
      pose.pose.orientation.y = q.y();
      pose.pose.orientation.z = q.z();
      pose.pose.orientation.w = q.w();

      path_msg.poses.push_back(pose);
    }

    // Publish the full predicted path
    predicted_path_pub_->publish(path_msg);

    //////////////////////////////////////////////////////////////
    // --- END OF: Publish the predicted trajectory for debugging ---
    //////////////////////////////////////////////////////////////

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

    //////////////////////////////////////////////////////////////
    // --- 7. Publish Prediction for Debugging ---
    //////////////////////////////////////////////////////////////
    // Extract the FIRST predicted state from the solution
    casadi::DM          p_next_dm  = solution.value(mpc_solver_.P_var_)(Slice(), 1); // col(1) is p[k+1]
    std::vector<double> p_next_vec = (std::vector<double>)p_next_dm;

    // Reconstruct the predicted nu[k+1] and eta[k+1] from p[k+1]
    Eigen::Map<Eigen::Vector<double, 9>> p_next_eigen(p_next_vec.data());
    Eigen::Vector<double, 6>             delta_x_next = p_next_eigen.head<6>();
    Eigen::Vector3d                      eta_next     = p_next_eigen.tail<3>();
    Eigen::Vector<double, 6>             x_next       = mpc_state_.x_k_ + delta_x_next;
    Eigen::Vector3d                      nu_next      = x_next.head<3>();

    // Create and populate the Odometry message
    nav_msgs::msg::Odometry prediction_msg;
    prediction_msg.header.stamp    = this->get_clock()->now();
    prediction_msg.header.frame_id = "map"; // Or your odom frame

    // Populate predicted pose (eta)
    prediction_msg.pose.pose.position.y = eta_next(0); // N
    prediction_msg.pose.pose.position.x = eta_next(1); // E
    tf2::Quaternion q;
    q.setRPY(0, 0, vk::wrapAngle(M_PI / 2.0 - eta_next(2))); // Convert back to ENU for plotting
    prediction_msg.pose.pose.orientation.x = q.x();
    prediction_msg.pose.pose.orientation.y = q.y();
    prediction_msg.pose.pose.orientation.z = q.z();
    prediction_msg.pose.pose.orientation.w = q.w();

    // Populate predicted twist (nu)
    prediction_msg.twist.twist.linear.x  = nu_next(0); // u
    prediction_msg.twist.twist.linear.y  = nu_next(1); // v
    prediction_msg.twist.twist.angular.z = nu_next(2); // r

    prediction_pub_->publish(prediction_msg);

    //////////////////////////////////////////////////////////////
    // --- END OF: 7. Publish Prediction for Debugging ---
    //////////////////////////////////////////////////////////////

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
    // casadi::MX cost = 0;
    casadi::MX cost = casadi::MX::zeros(1, 1);

    casadi::MX Q_mx = casadi::MX(mpc_params_.Q_);
    casadi::MX R_mx = casadi::MX(mpc_params_.R_);

    for (int i = 0; i < mpc_params_.P_; i++)
    {
      casadi::MX state_error       = mpc_solver_.P_var_(Slice(), i + 1) - mpc_solver_.p_ref_param_; // TODO wrap yaw
      casadi::MX yaw_error_raw     = state_error(8);
      casadi::MX yaw_error_wrapped = casadi::MX::atan2(casadi::MX::sin(yaw_error_raw), casadi::MX::cos(yaw_error_raw));
      state_error(8)               = yaw_error_wrapped;

      auto control_change = mpc_solver_.U_var_(Slice(), i);
      cost += casadi::MX::mtimes({state_error.T(), Q_mx, state_error});
      cost += casadi::MX::mtimes({control_change.T(), R_mx, control_change});
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
    // casadi::Dict opts;
    // opts["qpsol"]         = "qpoases";
    // opts["qpsol_options"] = casadi::Dict{{"printLevel", "none"}};
    // opts["print_time"]    = false;
    // opts["expand"]        = true;

    // opts["qpsol"]              = "qrqp"; // <-- change from "qpoases" to "osqp"
    // opts["qpsol_options"]      = casadi::Dict{{"verbose", false}};
    // opts["print_time"]         = true;
    // opts["sqpmethod.max_iter"] = 100;
    // opts["expand"]             = true;

    // casadi::Dict opts;
    // opts["qpsol"]              = "qrqp"; // or "osqp", "qpoases"
    // opts["qpsol_options"]      = casadi::Dict{{"verbose", false}};
    // opts["sqpmethod.max_iter"] = 100; // maximum iterations
    // opts["print_time"]         = true;
    // opts["expand"]             = true;

    // casadi::Function solver = nlpsol("solver", "sqpmethod", nlp, opts);

    // --- Create the Solver ---
    // Use IPOPT via the Opti interface and set max iterations
    casadi::Dict opts_ipopt;
    opts_ipopt["ipopt.max_iter"]    = 200; // <-- max iterations you want
    opts_ipopt["ipopt.print_level"] = 1;   // quiet IPOPT
    opts_ipopt["print_time"]        = false;
    opts_ipopt["expand"]            = false;

    // If you want a different QP solver for the inner QPs (used by some methods), you can set it like this:
    // opts_ipopt["qpsol"] = "qrqp";  // or "osqp", "qpoases" — optional

    // Attach IPOPT to the Opti instance
    mpc_solver_.opti_->solver("ipopt", opts_ipopt);
    RCLCPP_INFO(this->get_logger(), "MPC solver has been set up successfully with OSQP.");
  }
  void loadParameters()
  {
    MPCParams p_default;

    mpc_params_.P_  = this->declare_parameter<int>("P", p_default.P_);
    mpc_params_.Ts_ = this->declare_parameter<double>("Ts", p_default.Ts_);

    std::vector<double> Q_vec       = this->declare_parameter<std::vector<double>>("Q_weights");
    std::vector<double> R_vec       = this->declare_parameter<std::vector<double>>("R_weights");
    std::vector<double> tau_min_vec = this->declare_parameter<std::vector<double>>("tau_min");
    std::vector<double> tau_max_vec = this->declare_parameter<std::vector<double>>("tau_max");

    mpc_params_.Q_       = casadi::DM::diag(Q_vec);
    mpc_params_.R_       = casadi::DM::diag(R_vec);
    mpc_params_.tau_min_ = casadi::DM::reshape(tau_min_vec, 3, 1);
    mpc_params_.tau_max_ = casadi::DM::reshape(tau_max_vec, 3, 1);
  }
  // ROS2 Subscribers and Publishers
  rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr pos_cmd_sub_;
  rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr         state_sub_;
  rclcpp::Publisher<geometry_msgs::msg::Wrench>::SharedPtr         wrench_pub_;
  rclcpp::TimerBase::SharedPtr                                     timer_;
  rclcpp::Publisher<nav_msgs::msg::Path>::SharedPtr                predicted_path_pub_;
  rclcpp::Publisher<nav_msgs::msg::Odometry>::SharedPtr            prediction_pub_;

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
