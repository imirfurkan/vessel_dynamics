#pragma once
#include <cmath>
#include <Eigen/Dense>

namespace vk // TODO change name
{

inline constexpr double dt   = 0.02;     // 50 Hz
inline constexpr double m11_ = 2389.657; // [kg]
inline constexpr double m22_ = 2533.911; // [kg]
inline constexpr double m23_ = 62.386;   // [kg]
inline constexpr double m32_ = 28.141;   // [kg]
inline constexpr double m33_ = 5068.910; // [kg] or [kg*m2]

// ---------- Damping matrix D(nu) coefficients (fully coupled model) ----------
// Surge
inline constexpr double Xu_   = -27.632;  // [kg/s]
inline constexpr double Xuu_  = -110.064; // [kg/s]
inline constexpr double Xuuu_ = -13.965;  // [kg/s]
// Sway
inline constexpr double Yv_   = -52.947;   // [kg/s]
inline constexpr double Yvv_  = -116.486;  // [kg/s]
inline constexpr double Yvvv_ = -24.313;   // [kg/s]
inline constexpr double Yvr_  = 572.141;   // [kg/s]
inline constexpr double Yrv_  = -1540.383; // [kg/s]
inline constexpr double Yrr_  = -115.457;  // [kg/s]
inline constexpr double Yr_   = 24.372;    // [kg/s]
// Yaw
inline constexpr double Nv_   = 3.524;    // [kg/s]
inline constexpr double Nrv_  = 336.827;  // [kg/s]
inline constexpr double Nvv_  = -0.832;   // [kg/s]
inline constexpr double Nr_   = -122.860; // [kg/s]
inline constexpr double Nrr_  = -874.428; // [kg/s]
inline constexpr double Nvr_  = -121.957; // [kg/s]
inline constexpr double Nrrr_ = 0.0;      // [kg/s]

inline const Eigen::Matrix3d M_ = (Eigen::Matrix3d() <<

                                       m11_,
                                   0,
                                   0,

                                   0,
                                   m22_,
                                   m23_,

                                   0,
                                   m32_,
                                   m33_)
                                      .finished(); // Use .finished() for direct initialization with operator<<

inline double wrapAngle(double angle) // TODO learn inline
{
  constexpr double two_pi = 2.0 * M_PI;

  angle = std::fmod(angle, two_pi);

  if (angle < 0)
    angle += two_pi;

  if (angle > M_PI)
    angle -= two_pi;

  return angle;
}

inline Eigen::Vector3d clampVec3(const Eigen::Vector3d& x, const Eigen::Vector3d& lo, const Eigen::Vector3d& hi)
{
  return x.cwiseMax(lo).cwiseMin(hi);
}

inline Eigen::Matrix3d calculateRotationMatrix(const Eigen::Vector3d& pos)
{
  double cy = std::cos(pos(2)), sy = std::sin(pos(2)); // pos(2) is yaw.

  Eigen::Matrix3d Jnb;
  // clang-format off
  Jnb << cy, -sy, 0,
         sy,  cy, 0,
          0,   0, 1;
  // clang-format on

  return Jnb;
}

// velocity dependent → compute dynamically
inline void updateCoriolisMatrix(const Eigen::Vector3d& nu, Eigen::Matrix3d& C_out)
{
  double u = nu(0), v = nu(1), r = nu(2);

  double c13 = -m22_ * v - m23_ * r;
  double c23 = m11_ * u;

  // clang-format off
  C_out <<  0.0,  0.0, c13,
            0.0,  0.0, c23,
           -c13, -c23, 0.0;
  // clang-format on
}

// velocity dependent → compute dynamically
inline void updateDampingMatrix(const Eigen::Vector3d& nu, Eigen::Matrix3d& D_out)
{
  double u = nu(0), v = nu(1), r = nu(2);

  double d11 = -Xu_ - Xuu_ * std::abs(u) - Xuuu_ * u * u;
  double d22 = -Yv_ - Yvv_ * std::abs(v) - Yrv_ * std::abs(r) - Yvvv_ * v * v;
  double d23 = -Yr_ - Yvr_ * std::abs(v) - Yrr_ * std::abs(r);
  double d32 = -Nv_ - Nvv_ * std::abs(v) - Nrv_ * std::abs(r);
  double d33 = -Nr_ - Nvr_ * std::abs(v) - Nrr_ * std::abs(r) - Nrrr_ * r * r;

  // clang-format off
  D_out << d11, 0.0,  0.0,
           0.0, d22,  d23,
           0.0, d32,  d33;
  // clang-format on
}

} // namespace vk