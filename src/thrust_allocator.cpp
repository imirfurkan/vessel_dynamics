#include "vessel_kinematics/thrust_allocator.hpp"
#include <algorithm> // TODO when is this needed
#include <cmath>
#include <vector>
#include <chrono>
#include <osqp.h>
#include <Eigen/Dense> // TODO removing this still builds but it takes 4x time time
#include <Eigen/Sparse>
#include <iostream>
#include "vessel_kinematics/utils.hpp"

// Helper function declarations (at the bottom)
static inline Eigen::Vector2d clampVec(const Eigen::Vector2d& x, const Eigen::Vector2d& lo, const Eigen::Vector2d& hi);

Eigen::Matrix<double, 3, 2> Tmatrix(const Eigen::Vector2d& a, double Lx)
{
  const double c1 = std::cos(a(0)), s1 = std::sin(a(0));
  const double c2 = std::cos(a(1)), s2 = std::sin(a(1));

  Eigen::Matrix<double, 3, 2> T;

  // clang-format off
  T <<      c1,       c2,
            s1,       s2,
       Lx * s1, -Lx * s2;
  // clang-format on
  return T;
}

void dT_dalpha(const Eigen::Vector2d& a, double Lx, Eigen::Matrix<double, 3, 2>& dT1, Eigen::Matrix<double, 3, 2>& dT2)
{
  // Note: This function uses "output parameters" (dT1, dT2) passed by
  // reference. This is more efficient than returning large matrix objects,
  // as it avoids making unnecessary copies.
  const double c1 = std::cos(a(0)), s1 = std::sin(a(0));
  const double c2 = std::cos(a(1)), s2 = std::sin(a(1));

  // partial derivative of the T matrix wrt azimuth angle 1
  // clang-format off
  dT1 <<     -s1, 0.0,
              c1, 0.0,
         Lx * c1, 0.0;

  // partial derivative of the T matrix wrt azimuth angle 1
  dT2 << 0.0,    -s2,
         0.0,     c2,
         0.0, -Lx*c2;
  // clang-format on
}

// clang-format off
// --- Custom Deleters for OSQP C-API objects (for RAII) ---
struct OSQPSettingsDeleter { void operator()(OSQPSettings *s) const { if (s) c_free(s); } };
struct OSQPDataDeleter { void operator()(OSQPData *d) const { if(d) { if (d->P) c_free(d->P); if (d->A) c_free(d->A); c_free(d); } } };
struct OSQPWorkspaceDeleter { void operator()(OSQPWorkspace *w) const { if (w) osqp_cleanup(w); } };
// clang-format on

// main dispatcher function
TAResult allocate_tau(const Eigen::Vector3d& tau_des, const TAState& state, const TAParams& P, double dt)
{
  // --- QP Problem Definition ---
  // Solves for x = [f_0, f_1, s_x, s_y, s_n, dα_0, dα_1]
  constexpr int NUM_VARS        = 7; // 2 thrusts, 3 slacks, 2 angle changes
  constexpr int NUM_CONSTRAINTS = 9;

  // --- Hessian Matrix P (Quadratic Cost) ---
  Eigen::VectorXd P_diag(NUM_VARS); // diagonal parts for the P matrix
  P_diag << P.Wf, P.Qs, P.Oa;
  P_diag *= 2.0; // OSQP minimizes 0.5 x'Px + q'x

  // Pre-allocate memory for efficiency. This avoids slow resizing later.
  std::vector<c_float> P_data;
  P_data.reserve(NUM_VARS);
  std::vector<c_int> P_i;
  P_i.reserve(NUM_VARS);
  std::vector<c_int> P_p;
  P_p.reserve(NUM_VARS + 1);

  // This loop fills the vectors for our diagonal matrix.
  for (int k = 0; k < NUM_VARS; ++k)
  {
    P_p.push_back(k);            // The data for column 'k' starts at index 'k'.
    P_i.push_back(k);            // The non-zero value in column 'k' is also in row 'k'.
    P_data.push_back(P_diag(k)); // The actual value is the k-th weight from our P_diag vector.
  }
  // The CSC format requires one final entry in P_p to mark the end.
  P_p.push_back(NUM_VARS);

  // --- Linear Cost q (zero) ---
  std::vector<c_float> q(NUM_VARS, 0.0); // we don't have a linear cost, we just minimize 0.5 x'Px

  const Eigen::Matrix<double, 3, 2> T0 = Tmatrix(state.alpha, P.Lx);

  Eigen::Matrix<double, 3, 2> dT1, dT2;
  dT_dalpha(state.alpha, P.Lx, dT1, dT2);

  // Linearization term for dα_0 and dα_1.
  Eigen::Vector3d       epsilon_     = Eigen::Vector3d::Constant(0.01); // unstick the system
  const Eigen::Vector3d term_dalpha0 = dT1.col(0) * state.f(0) + epsilon_;
  const Eigen::Vector3d term_dalpha1 = dT2.col(1) * state.f(1) + epsilon_;

  // Using a vector of triplets to build the sparse matrix A in a clean way.
  // We'll convert this to CSC format later.
  // A triplet is a simple (row, column, value) tuple.
  std::vector<Eigen::Triplet<double>> triplets; //  OSQP library is written in C and uses the more memory-efficient CSC
                                                //  format, which is not directly compatible with Eigen's triplet
                                                //  representation so we have to make it vector.

  triplets.reserve(21); // TODO what is this

  // Row 0-2: Physics constraint (T0*f + (dT/dalpha)*f_prev*dalpha + s = tau_des)
  // [f0]
  triplets.emplace_back(0, 0, T0(0, 0));
  triplets.emplace_back(1, 0, T0(1, 0));
  triplets.emplace_back(2, 0, T0(2, 0));
  // [f1]
  triplets.emplace_back(0, 1, T0(0, 1));
  triplets.emplace_back(1, 1, T0(1, 1));
  triplets.emplace_back(2, 1, T0(2, 1));
  // [sx, sy, sn]
  triplets.emplace_back(0, 2, -1.0);
  triplets.emplace_back(1, 3, -1.0);
  triplets.emplace_back(2, 4, -1.0);
  // [dα0]
  triplets.emplace_back(0, 5, term_dalpha0(0));
  triplets.emplace_back(1, 5, term_dalpha0(1));
  triplets.emplace_back(2, 5, term_dalpha0(2));
  // [dα1]
  triplets.emplace_back(0, 6, term_dalpha1(0));
  triplets.emplace_back(1, 6, term_dalpha1(1));
  triplets.emplace_back(2, 6, term_dalpha1(2));
  // Row 3: f0 lower bound
  triplets.emplace_back(3, 0, 1.0);
  // Row 4: f1 lower bound
  triplets.emplace_back(4, 1, 1.0);
  // Row 5: dα0 lower bound
  triplets.emplace_back(5, 5, 1.0);
  // Row 6: dα1 lower bound
  triplets.emplace_back(6, 6, 1.0);

  // NEW: Row 7-8: Force rate-of-change constraints
  // This constrains the final force f, not the change df.
  // The constraint is: f_prev - max_rate*dt <= f_current <= f_prev + max_rate*dt
  triplets.emplace_back(7, 0, 1.0); // for f0
  triplets.emplace_back(8, 1, 1.0); // for f1

  // Convert triplets to OSQP's CSC format. This is the most reliable way.
  Eigen::SparseMatrix<double> A_sparse(NUM_CONSTRAINTS, NUM_VARS);
  A_sparse.setFromTriplets(triplets.begin(), triplets.end());
  A_sparse.makeCompressed();

  // Extract CSC data
  std::vector<c_float> A_data(A_sparse.valuePtr(), A_sparse.valuePtr() + A_sparse.nonZeros());
  std::vector<c_int>   A_i(A_sparse.innerIndexPtr(), A_sparse.innerIndexPtr() + A_sparse.nonZeros());
  std::vector<c_int>   A_p(A_sparse.outerIndexPtr(), A_sparse.outerIndexPtr() + (A_sparse.outerSize() + 1));

  // --- Constraint Bounds l and u ---
  std::vector<c_float>  l(NUM_CONSTRAINTS), u(NUM_CONSTRAINTS);
  const Eigen::Vector2d dalpha_min       = P.dalpha_min_rate_rad_s * dt;
  const Eigen::Vector2d dalpha_max       = P.dalpha_max_rate_rad_s * dt;
  const double          max_force_change = P.max_force_rate * dt;

  // Physics constraint: T0*f + s = tau_des
  l[0] = u[0] = tau_des(0);
  l[1] = u[1] = tau_des(1);
  l[2] = u[2] = tau_des(2);
  // Thrust bounds
  l[3] = P.f_min(0), u[3] = P.f_max(0);
  l[4] = P.f_min(1), u[4] = P.f_max(1);
  // Asymmetric angle rate limits
  l[5] = dalpha_min(0), u[5] = dalpha_max(0);
  l[6] = dalpha_min(1), u[6] = dalpha_max(1);
  // NEW: Row 7-8: Force rate-of-change bounds
  l[7] = state.f(0) - max_force_change;
  u[7] = state.f(0) + max_force_change;
  l[8] = state.f(1) - max_force_change;
  u[8] = state.f(1) + max_force_change;

  // --- Setup and Solve with OSQP (using RAII for safety) --- // TODO detayli incelemedim
  std::unique_ptr<OSQPSettings, OSQPSettingsDeleter> settings((OSQPSettings*)c_malloc(sizeof(OSQPSettings)));
  osqp_set_default_settings(settings.get());
  settings->verbose  = 0;
  settings->polish   = true;
  settings->eps_abs  = 1e-5; // Tighter absolute tolerance
  settings->eps_rel  = 1e-5;
  settings->max_iter = 8000;

  std::unique_ptr<OSQPData, OSQPDataDeleter> data((OSQPData*)c_malloc(sizeof(OSQPData)));
  data->n = NUM_VARS;
  data->m = NUM_CONSTRAINTS;
  data->P = csc_matrix(data->n, data->n, P_data.size(), P_data.data(), P_i.data(), P_p.data());
  data->q = q.data();
  data->A = csc_matrix(data->m, data->n, A_data.size(), A_data.data(), A_i.data(), A_p.data());
  data->l = l.data();
  data->u = u.data();

  OSQPWorkspace* raw_workspace = nullptr;
  osqp_setup(&raw_workspace, data.get(), settings.get());
  std::unique_ptr<OSQPWorkspace, OSQPWorkspaceDeleter> work(raw_workspace);

  // Record the start time using std::chrono::steady_clock
  auto start_time = std::chrono::steady_clock::now();

  // The line you want to time
  osqp_solve(work.get());

  // Record the end time
  auto end_time = std::chrono::steady_clock::now();

  // Calculate the duration
  std::chrono::duration<double> duration = end_time - start_time;

  // Print the elapsed time
  std::cout << "OSQP solver took: " << duration.count() << " seconds." << std::endl;

  // --- Extract Results ---
  TAResult   out;
  const bool is_solved = work && work->info->status_val == OSQP_SOLVED;
  out.success          = is_solved;
  if (is_solved)
  {
    const c_float*  sol_x = work->solution->x;
    Eigen::Vector2d f_sol(sol_x[0], sol_x[1]);
    Eigen::Vector2d da_sol(sol_x[5], sol_x[6]);

    out.alpha    = state.alpha + da_sol;
    out.alpha(0) = vk::wrapAngle(out.alpha(0));
    out.alpha(1) = vk::wrapAngle(out.alpha(1));
    // out.f        = clampVec(f_sol, P.f_min, P.f_max);
    out.f = f_sol; // TODO let's not clamp it and trust the solver
    std::cout << "azimuth 0 change:" << da_sol[0] * 180 / M_PI << std::endl;
  }
  else
  {
    // If solver fails, return a safe state (no change in angle, zero thrust).
    // The caller MUST check out success to handle this failure case.
    out.alpha = state.alpha;
    out.f.setZero();
    // out.f = state.f; // lets keep it at where it was
    std::cout << "solver failed" << std::endl;
  }
  return out;
}

static inline Eigen::Vector2d clampVec(const Eigen::Vector2d& x, const Eigen::Vector2d& lo, const Eigen::Vector2d& hi)
{
  return x.cwiseMax(lo).cwiseMin(hi);
}