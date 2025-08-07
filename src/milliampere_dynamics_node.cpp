#include <memory>
#include <cmath>
#include <rclcpp/rclcpp.hpp>

#include <geometry_msgs/msg/wrench.hpp>
#include <nav_msgs/msg/odometry.hpp>
#include <geometry_msgs/msg/transform_stamped.hpp>
#include <tf2_ros/transform_broadcaster.h>
#include <tf2/LinearMath/Quaternion.h>

#include <Eigen/Dense>

using namespace std::chrono_literals;

class MilliampereDynamics : public rclcpp::Node
{
public:
  MilliampereDynamics()              // constructor of the class
      : Node("milliampere_dynamics") // initializer list that calls the base class class constructor
                                     // (rclcpp::Node) with the name "milliampere_dynamics"/
  {
    //----------------------------------------------------------------
    // 1) Parameters
    //----------------------------------------------------------------
    // clang-format off
    M_ << m11_,    0,    0,
             0, m22_, m23_,
             0, m32_, m33_;
    // clang-format on

    //----------------------------------------------------------------
    // 2) State initialization
    //----------------------------------------------------------------
    nu_.setZero();  // [u, v, r]
    eta_.setZero(); // [x, y, ψ]
    tau_.setZero(); // [X, Y, N]
  }

private:
  void step()
  {
    const double dt = 0.02; // 50 Hz

    // ------------------------------------------------------------
    //  A) Rigid‐body dynamics
    // ------------------------------------------------------------
    double ca1 = std::cos(azimuthAngles_(0)), ca2 = std::cos(azimuthAngles_(1));
    double sa1 = std::sin(azimuthAngles_(0)), sa2 = std::sin(azimuthAngles_(1));

    // clang-format off
    T_ <<       ca1,      ca2,
                sa1,      sa2,
          Lx_ * sa1, Lx * sa2;
    // clang-format on

    tau_ = T_ * thrustforces_;

    // Mν̇ + C(v)v + D(v)v = τ
    Eigen::Vector3d nu_dot = (-C_ * nu_ - D_ * nu_ + tau_) / M_;

    // integrate body‐frame velocity νk+1​=νk​+ν˙k​⋅Δt
    nu_ += nu_dot * dt;

    // ------------------------------------------------------------
    //  B) Kinematic mapping to inertial pose
    // ------------------------------------------------------------

    double cy = std::cos(eta_(2)), sy = std::sin(eta_(2));

    // clang-format off
    Jnb << cy, -sy, 0,
           sy,  cy, 0,
            0,   0, 1;
    // clang-format on

    Eigen::Vector3d eta_dot = Jnb * nu_;
    eta_ += eta_dot * dt;
  }

  // ---------- node state & parameters ------------
  double          m_  = 1800;     // mass [kg]
  double          Lx_ = 2;        // distance from CO to thrusters [m] TODO find the actual value
  Eigen::Matrix3d Jnb;            // rotation matrix from body to inertial
  Eigen::Matrix3d M_;             // mass inertia matrix
  Eigen::Matrix3d C_;             // coriolis-centripetal matrix
  Eigen::Matrix3d D_;             // damping matrix
  Eigen::Vector3d eta_;           // position [x, y, ψ]
  Eigen::Vector3d nu_;            // body-frame twist [u, v, r]
  Eigen::Vector3d tau_;           // [X, Y, N]
  Eigen::Vector2d azimuthAngles_; // [α1, α2]
  Eigen::Vector2d thrustForces_;  // [F1, F2]
  Eigen::Matrix<double, 3, 2> T_; // Thrust configuration matrix

  double m11_ = 2389.657; // [kg]
  double m22_ = 2533.911; // [kg]
  double m23_ = 62.386;   // [kg]
  double m32_ = 28.141;   // [kg]
  double m33_ = 5068.910; // [kg]

  // ---------- Damping matrix D(nu) coefficients (fully coupled model) ----------
  // Surge
  double Xu_   = -27.632;  // [kg/s]
  double Xuu_  = -110.064; // [kg/s]
  double Xuuu_ = -13.965;  // [kg/s]
  // Sway
  double Yv_   = -52.947;   // [kg/s]
  double Yvv_  = -116.486;  // [kg/s]
  double Yvvv_ = -24.313;   // [kg/s]
  double Yvr_  = 572.141;   // [kg/s]
  double Yrv_  = -1540.383; // [kg/s]
  double Yrr_  = -115.457;  // [kg/s]
  double Yr_   = 24.372;    // [kg/s]
  // Yaw
  double Nv_   = 3.524;    // [kg/s]
  double Nrv_  = 336.827;  // [kg/s]
  double Nvv_  = -0.832;   // [kg/s]
  double Nr_   = -122.860; // [kg/s]
  double Nrr_  = -874.428; // [kg/s]
  double Nvr_  = -121.957; // [kg/s]
  double Nrrr_ = 0.0;      // [kg/s]

  // Coriolis matrix is velocity dependent → compute dynamically
  void updateCoriolisMatrix(const Eigen::Vector3d& nu)
  {
    double u = nu(0), v = nu(1), r = nu(2);

    double c13 = -m22_ * v - m23_ * r;
    double c23 = m11_ * u;

    // clang-format off
    C_ <<  0.0,  0.0, c13,
           0.0,  0.0, c23,
          -c13, -c23, 0.0;
    // clang-format on
  }

  void updateDampingMatrix(const Eigen::Vector3d& nu)
  {
    double u = nu(0), v = nu(1), r = nu(2);

    double d11 = -Xu_ - Xuu_ * std::abs(u) - Xuuu_ * u * u;
    double d22 = -Yv_ - Yvv_ * std::abs(v) - Yrv_ * std::abs(r) - Yvvv_ * v * v;
    double d23 = -Yr_ - Yvr_ * std::abs(v) - Yrr_ * std::abs(r);
    double d32 = -Nv_ - Nvv_ * std::abs(v) - Nrv_ * std::abs(r);
    double d33 = -Nr_ - Nvr_ * std::abs(v) - Nrr_ * std::abs(r) - Nrrr_ * r * r;

    // clang-format off
    D_ << d11, 0.0,  0.0,
          0.0, d22,  d23,
          0.0, d32,  d33;
    // clang-format on
  }
}