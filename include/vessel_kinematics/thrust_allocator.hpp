#pragma once

#include <Eigen/Dense>
#include <memory> // Required for std::unique_ptr // TODO check that function

// clang-format off

// --- Forward declarations for OSQP C-style structs --- // TODO check
struct OSQPSettings;
struct OSQPWorkspace;
struct OSQPData;

/**
 * @file thrust_allocator.hpp
 * @brief An optimal thruster allocator for a 2-thruster, 3-DOF vessel.
 *
 * This allocator computes individual thruster angles (alpha) and forces (f)
 * to achieve a desired generalized force vector tau = [Fx, Fy, N].
 * It uses a Quadratic Programming (QP) approach for optimal results.
 */

struct TAParams
{
  double Lx = 1.8; // [m]
  const double dt = 0.02; // 50 Hz

  // TODO why this initialization method
  Eigen::Vector2d f_min{0.0, 0.0};     // minimum thrust force per thruster [N]
  Eigen::Vector2d f_max{500.0, 500.0}; // maximum thrust force per thruster [N]

  // Azimuth rate limits ~ From Table 3.9 of Pedersen
  // Index 0 is for thruster 1, index 1 is for thruster 2.
  Eigen::Vector2d dalpha_min_rate_rad_s{-34.46, -37.53};
  Eigen::Vector2d dalpha_max_rate_rad_s{34.46, 37.53};

  // --- QP Optimization Weights ---
  // These weights define the "priorities" of the optimization.
  // Tweak these to change the allocator's behavior.
  // These vectors are later concatenated into one large diagonal vector
  // to build the sparse matrix 'P' that OSQP needs, inside osqp_path.
  Eigen::Vector2d Wf{1.0, 1.0};      // Higher values prioritize using less thrust.
  Eigen::Vector3d Qs{1e3, 1e3, 1e3}; // This should be high to ensure the desired tau is met.
  Eigen::Vector2d Oa{10.0, 10.0} // Higher values result in smoother, less aggressive angle changes.
};

/// @brief Current state of the thruster system. This is updated at each timestep.
struct TAState
{
  Eigen::Vector2d alpha = Eigen::Vector2d::Zero(); // Current thruster angles [rad, rad].
                                                   // TODO how is this updated each time step
};

/// @brief The result computed by the thruster allocator.
struct TAResult
{
  Eigen::Vector2d alpha;
  Eigen::Vector2d f;
  bool success = false; // True if the OSQP solver found an optimal solution. If false, the returned command is a safe state (current angles and zero thrust).
};

Eigen::Matrix<double, 3, 2> Tmatrix(const Eigen::Vector2d& a, double Lx);

void dT_dalpha(const Eigen::Vector2d&       a,
               double                       Lx,
               Eigen::Matrix<double, 3, 2>& dT1,
               Eigen::Matrix<double, 3, 2>& dT2);


/// @brief Main entry point to compute thruster commands using OSQP.
/// @return A TAResult struct. The 'success' flag MUST be checked by the caller. // TODO is that so
TAResult allocate_tau(const Eigen::Vector3d &tau_des,
                      const TAState &state,
                      const TAParams &P,
                      double dt);

