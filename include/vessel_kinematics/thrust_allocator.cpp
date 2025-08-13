#include "thrust_allocator.hpp"
#include <algorithm> // TODO when is this needed
#include <cmath>
#include <vector>
#include <osqp.h>

// Helper function declarations (at the bottom)
static inline Vec2   clampVec(const Vec2& x, const Vec2& lo, const Vec2& hi);
static inline double wrapAngle(double angle);

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
}

void dT_dalpha(const Eigen::Vector2d&       a,
               double                       Lx,
               Eigen::Matrix<double, 3, 2>& dT1,
               Eigen::Matrix<double, 3, 2>& dT2)
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
TAResult
allocate_tau(const Eigen::Vector3d& tau_des, const TAState& state, const TAParams& P, double dt)
{
  // --- QP Problem Definition ---
  // Solves for x = [f_0, f_1, s_x, s_y, s_n, dα_0, dα_1]
  constexpr int NUM_VARS        = 7; // 2 thrusts, 3 slacks, 2 angle changes
  constexpr int NUM_CONSTRAINTS = 7;

  // --- Hessian Matrix P (Quadratic Cost) ---
  Eigen::VectorXd P_diag(NUM_VARS); // diagonal parts for the P matrix
  P_diag << P.Wf, P.Qs, P.Oa;
  P_diag *= 2.0; // OSQP minimizes 0.5 x'Px + q'x

  // Pre-allocate memory for efficiency. This avoids slow resizing later.
  P_data.reserve(NUM_VARS);
  P_i.reserve(NUM_VARS);
  P_p.reserve(NUM_VARS + 1); // accomply with CSC data structure

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

  const Eigen::Matrix<double, 3, 2> T0     = Tmatrix(state.alpha, P.Lx);
  std::vector<c_float>              A_data = {
      T0(0, 0),
      T0(1, 0),
      T0(2, 0),
      1.0,
      T0(0, 1),
      T0(1, 1),
      T0(2, 1),
      1.0,
      -1.0,
      -1.0,
      -1.0,
      1.0,
      1.0}; // list of all the non-zero values from the A matrix, read column by column.
  std::vector<c_int> A_i = {0, 1, 2, 3, 0, 1, 2, 4, 0, 1, 2, 5, 6};
  std::vector<c_int> A_p = {0, 4, 8, 9, 10, 11, 12, 13};

  // --- Constraint Bounds l and u ---
  std::vector<c_float>  l(NUM_CONSTRAINTS), u(NUM_CONSTRAINTS);
  const Eigen::Vector2d dalpha_min = P.dalpha_min_rate_rad_s * P.dt;
  const Eigen::Vector2d dalpha_max = P.dalpha_max_rate_rad_s * P.dt;

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

  // --- Setup and Solve with OSQP (using RAII for safety) --- // TODO detayli incelemedim
  std::unique_ptr<OSQPSettings, OSQPSettingsDeleter> settings(
      (OSQPSettings*)c_malloc(sizeof(OSQPSettings)));
  osqp_set_default_settings(settings.get());
  settings->verbose = 0;

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

  osqp_solve(work.get());

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
    out.alpha(0) = wrapAngle(out.alpha(0));
    out.alpha(1) = wrapAngle(out.alpha(1));
    out.f        = clampVec(f_sol, P.f_min, P.f_max);
  }
  else
  {
    // If solver fails, return a safe state (no change in angle, zero thrust).
    // The caller MUST check out success to handle this failure case.
    out.alpha = state.alpha;
    out.f.setZero();
  }
  return out;
}

static inline Vec2 clampVec(const Vec2& x, const Vec2& lo, const Vec2& hi)
{
  return x.cwiseMax(lo).cwiseMin(hi);
}

static inline double wrapAngle(double angle)
{
  const double two_pi = 2.0 * M_PI;

  // Step 1: First, wrap the angle into the [0, 2π) range.
  // After this block, 'angle' is guaranteed to be positive.
  angle = std::fmod(angle, two_pi);
  if (angle < 0)
    angle += two_pi;

  // Step 2: Now, move angles from the top half (π, 2π) to the equivalent
  // negative half (-π, 0).
  if (angle > M_PI)
    angle -= two_pi;

  return angle;
}