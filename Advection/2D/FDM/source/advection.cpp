#include <advection.h>
#include <utility.h>

#include <iostream>
#include <cmath>

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void upwindMethod (Ref<VectorXd> u,
                   const int N,
                   const Vector2d v,
                   const double delta_t) {
  const double h     = 1.0 / (N + 1),
               sigma_x = v[0] * delta_t / h,
               sigma_y = v[1] * delta_t / h;

  VectorXd u1 = u;
  for (int i = 0; i < N; ++i) 
    for (int j = 0; j < N; ++j)
      u[i * N + j] = bc(u1, N, i, j) + sigma_x * (bc(u1, N, i - 1, j) - bc (u1, N, i, j)) + sigma_y * (bc (u1, N, i, j - 1) - bc (u1, N, i, j));
}

void frommMethod (Ref<VectorXd> u,
                  const int N,
                  const Vector2d v,
                  const double delta_t) {
  const double h       = 1.0 / (N + 1),
               sigma_x = v[0] * delta_t / h,
               sigma_y = v[1] * delta_t / h;

  VectorXd u1 = u;
  for (int i = 0; i < N; ++i) 
    for (int j = 0; j < N; ++j) 
      u[i * N + j] = bc (u1, N, i, j) + sigma_x * 
                           ((bc (u1, N, i - 1, j) + (1.0 - sigma_x) * 0.25 * (bc (u1, N, i, j) - bc (u1, N, i - 2, j))) - 
                            (bc (u1, N, i, j) + (1.0 - sigma_x) * 0.25 * (bc (u1, N, i + 1, j) - bc (u1, N, i - 1, j)))) +
                                        sigma_y *
                           ((bc (u1, N, i, j - 1) + (1.0 - sigma_y) * 0.25 * (bc (u1, N, i, j) - bc (u1, N, i, j - 2))) -
                            (bc (u1, N, i, j) + (1.0 - sigma_y) * 0.25 * (bc (u1, N, i, j + 1) - bc (u1, N, i, j - 1))));
}

void frommVanLeer (Ref<VectorXd> u,
                   const int N,
                   const Vector2d v,
                   const double delta_t) {
  const double h       = 1.0 / (N + 1),
               sigma_x = v[0] * delta_t / h,
               sigma_y = v[1] * delta_t / h;

  VectorXd u1 = u;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      double thetaL_x = min3 (2.0 * std::abs (bc (u1, N, i - 1, j) - bc (u1, N, i - 2, j)),
                              0.5 * std::abs (bc (u1, N, i, j) - bc (u1, N, i - 2, j)),
                              2.0 * std::abs (bc (u1, N, i, j) - bc (u1, N, i - 1, j)));
      double thetaR_x = min3 (2.0 * std::abs (bc (u1, N, i, j) - bc (u1, N, i - 1, j)),
                              0.5 * std::abs (bc (u1, N, i + 1, j) - bc (u1, N, i - 1, j)),
                              2.0 * std::abs (bc (u1, N, i + 1, j) - bc (u1, N, i, j)));
      double thetaL_y = min3 (2.0 * std::abs (bc (u1, N, i, j - 1) - bc (u1, N, i, j - 2)),
                              0.5 * std::abs (bc (u1, N, i, j) - bc (u1, N, i, j - 2)),
                              2.0 * std::abs (bc (u1, N, i, j) - bc (u1, N, i, j - 1)));
      double thetaR_y = min3 (2.0 * std::abs (bc (u1, N, i, j) - bc (u1, N, i, j - 1)),
                              0.5 * std::abs (bc (u1, N, i, j + 1) - bc (u1, N, i, j - 1)),
                              2.0 * std::abs (bc (u1, N, i, j + 1) - bc (u1, N, i, j)));
      double phiL_x = (bc (u1, N, i, j) - bc (u1, N, i - 1, j)) *
                      (bc (u1, N, i - 1, j) - bc (u1, N, i - 2, j));
      double phiR_x = (bc (u1, N, i + 1, j) - bc (u1, N, i, j)) *
                      (bc (u1, N, i, j) - bc (u1, N, i - 1, j));
      double phiL_y = (bc (u1, N, i, j) - bc (u1, N, i, j - 1)) *
                      (bc (u1, N, i, j - 1) - bc (u1, N, i, j - 2));
      double phiR_y = (bc (u1, N, i, j + 1) - bc (u1, N, i, j)) *
                      (bc (u1, N, i, j) - bc (u1, N, i, j - 1));
      double VLDeltaL_x = (phiL_x > 0) ? 
                            copysign (thetaL_x, 
                                      (bc (u1, N, i, j) - bc (u1, N, i - 2, j)))
                            : 0;
      double VLDeltaR_x = (phiR_x > 0) ?
                            copysign (thetaR_x,
                                      (bc (u1, N, i + 1, j) - bc (u1, N, i - 1, j)))
                            : 0;
      double VLDeltaL_y = (phiL_y > 0) ?
                            copysign (thetaL_y,
                                      (bc (u1, N, i, j) - bc (u1, N, i, j - 2)))
                            : 0;
      double VLDeltaR_y = (phiR_y > 0) ?
                            copysign (thetaR_y,
                                      (bc (u1, N, i, j + 1) - bc (u1, N, i, j - 1)))
                            : 0;
      double fluxL_x = bc (u1, N, i - 1, j) + 0.5 * (1.0 - sigma_x) * VLDeltaL_x;
      double fluxR_x = bc (u1, N, i, j) + 0.5 * (1.0 - sigma_x) * VLDeltaR_x;
      double fluxL_y = bc (u1, N, i, j - 1) + 0.5 * (1.0 - sigma_y) * VLDeltaL_y;
      double fluxR_y = bc (u1, N, i, j) + 0.5 * (1.0 - sigma_y) * VLDeltaR_y;
      u[i * N + j] = bc (u1, N, i, j) + sigma_x * (fluxL_x - fluxR_x) 
                                      + sigma_y * (fluxL_y - fluxR_y);
    }
  }
}

