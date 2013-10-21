#include <advection.h>
#include <utility.h>

#include <iostream>
#include <cmath>

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void upwindMethod (Ref<VectorXd> u,
                   const double dt,
                   const double h,
                   const VectorXd v) {
  const int N     = u.rows ();
  double    sigmaL;
  double    sigmaR;

  VectorXd u1 = u;
  for (int i = 0; i < int (N); ++i) {
    sigmaR = sigmaL = 0;
    u[i] = bc(u1, N, i) - 1 / h * (v[0] * dt * (bc(u1, N, i) - bc(u1, N, i-1)) + (v[0] * h / 2 * h - v[0] * v[0] * dt * dt / 2) * (sigmaR - sigmaL));
  }
}

void frommMethod (Ref<VectorXd> u,
                  const double delta_t,
                  const double h,
                  const VectorXd v) {
  const int    N     = u.rows ();
  const double sigma = v[0] * delta_t / h;

  VectorXd u1 = u;
  for (int i = 0; i < int (N); ++i) 
    u[i] = u1[i] + sigma *
                   ((bc (u1, N, i - 1) + (1.0 - sigma) * 0.25 * (bc (u1, N, i) - bc (u1, N, i - 2))) -
                    (bc (u1, N, i) + (1.0 - sigma) * 0.25 * (bc (u1, N, i + 1) - bc (u1, N, i - 1))));
}


void frommVanLeer (Ref<VectorXd> u,
                   const double delta_t,
                   const double h,
                   const VectorXd v) {
  double thetaL, thetaR, phiL, phiR, VLDeltaL, VLDeltaR, fluxL, fluxR;
  
  const int    N     = u.rows ();
  const double sigma = v[0] * delta_t / h;

  VectorXd u1 = u;
  for (int i = 0; i < N; ++i) {
    thetaL = min3 (2.0 * std::abs (bc (u1, N, i - 1) - bc (u1, N, i - 2)),
                   0.5 * std::abs (u1[i] - bc (u1, N, i - 2)),
                   2.0 * std::abs (u1[i] - bc (u1, N, i - 1)) );
    thetaR = min3 (2.0 * std::abs (u1[i] - bc (u1, N, i - 1)),
                   0.5 * std::abs (bc (u1, N, i + 1) - bc (u1, N, i - 1)),
                   2.0 * std::abs (bc (u1, N, i + 1) - u1[i]) );
    phiL = (u1[i] - bc (u1, N, i - 1)) *
           (bc (u1, N, i - 1) - bc (u1, N, i - 2));
    phiR = (bc (u1, N, i + 1) - u1[i]) *
           (u1[i] - bc (u1, N, i - 1));
    VLDeltaL = (phiL > 0) ? copysign (thetaL, (u1[i] - bc (u1, N, i - 2)) ) : 0;
    VLDeltaR = (phiR > 0) ? copysign (thetaR, (bc (u1, N, i + 1) - bc (u1, N, i - 1)) ) : 0;
    fluxL = bc (u1, N, i - 1) + 0.5 * (1.0 - sigma) * VLDeltaL;
    fluxR = u1[i] + 0.5 * (1.0 - sigma) * VLDeltaR;
    u[i] = u1[i] + sigma * (fluxL - fluxR);
  }
}
