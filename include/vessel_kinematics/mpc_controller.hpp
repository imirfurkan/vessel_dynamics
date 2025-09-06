#pragma once

#include <casadi/casadi.hpp>
#include <Eigen/Dense>

using namespace casadi;

// --- MPC Tuning Parameters ---
struct MPCParams
{
  int        P_;       // Prediction horizon (number of steps)
  double     Ts_;      // Control timestep (seconds), e.g., 0.1s for 10Hz
  casadi::DM Q_;       // Weighting matrix for state error (augmented)
  casadi::DM R_;       // Weighting matrix for control input changes
  casadi::DM tau_min_; // Minimum thruster force/torque
  casadi::DM tau_max_; // Maximum thruster force/torque
};

struct MPCState
{
  // --- System State & Model ---
  Eigen::Vector<double, 6> x_k_;           // Current full state [nu, eta]
  Eigen::Vector<double, 6> x_k_minus_1_;   // Previous full state (for delta_x)
  Eigen::Vector3d          tau_k_minus_1_; // Last applied control input (for delta_tau)
};

struct MPCSolver
{
  // --- CasADi Solver Objects ---
  std::shared_ptr<casadi::Opti> opti_; // The optimization problem stack

  // Symbolic VARIABLES the solver will find
  casadi::MX U_var_; // The sequence of control changes {Δτ_k, ..., Δτ_{k+P-1}}
  casadi::MX P_var_; // The sequence of predicted states {p_{k+1}, ..., p_{k+P}}

  // Symbolic PARAMETERS we will provide at runtime
  casadi::MX p0_param_;       // The initial augmented state p[k]
  casadi::MX p_ref_param_;    // The reference for the output y (part of p)
  casadi::MX tau_prev_param_; // The previous control input τ[k-1]

  // Matrices for the augmented model p[k+1] = A_bar*p[k] + B_bar*delta_tau[k]
  // These will be recalculated in every loop.
  casadi::MX A_bar_;
  casadi::MX B_bar_;
};