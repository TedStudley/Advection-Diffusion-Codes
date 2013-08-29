#include <advection.h>
#include <utility.h>

#include <iostream>
#include <cmath>

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void upwindMethod (Ref<VectorXd> u,
                   unsigned int N,
                   double v,
                   double delta_t) {
  const double h     = 1.0 / (N + 1),
               sigma = v * delta_t / h;

  if (sigma > 1) 
    cerr << "Warning! CFL condition " << sigma << " > 1. Numerical stability is not guaranteed." << endl;

  VectorXd u1 = u;
  for (int i = 0; i < int (N); ++i) { 
    u[i] = bc(u1, N, i) + sigma * ( bc(u1, N, i - 1) - bc(u1, N, i));
  }
}

void frommMethod (Ref<VectorXd> u,
                  unsigned int N,
                  double v,
                  double delta_t) {
  const double h     = 1.0 / (N + 1);
        double sigma = v * delta_t / h;

  if (sigma > 1)
    cerr << "Warning! CFL condition " << sigma << " > 1. Numerical stability is not guaranteed." << endl;

  VectorXd u1 = u;
  for (int i = 0; i < int (N); ++i) 
    u[i] = u1[i] + sigma *
                   ((bc (u1, N, i - 1) + (1.0 - sigma) * 0.25 * (bc (u1, N, i) - bc (u1, N, i - 2))) -
                    (bc (u1, N, i) + (1.0 - sigma) * 0.25 * (bc (u1, N, i + 1) - bc (u1, N, i - 1))));
}

void frommVanLeer (Ref<VectorXd> u,
                   unsigned int N,
                   double v,
                   double delta_t) {
  const double h     = 1.0 / (N + 1),
               sigma = v * delta_t / h;

  if (sigma > 1)
    cerr << "Warning! CFL condition " << sigma << " > 1. Numerical stability is not guaranteed." << endl;

  double thetaL, thetaR, phiL, phiR, VLDeltaL, VLDeltaR, fluxL, fluxR;
  
  VectorXd u1 = u;
  for (int i = 0; i < int(N); ++i) {
    thetaL = min3 (2.0 * std::abs (bc (u1, N, i - 1) - bc (u1, N, i - 2)),
                   0.5 * std::abs (bc (u1, N, i) - bc (u1, N, i - 2)),
                   2.0 * std::abs (bc (u1, N, i) - bc (u1, N, i - 1)));
    thetaR = min3 (2.0 * std::abs (bc (u1, N, i) - bc (u1, N, i - 1)),
                   0.5 * std::abs (bc (u1, N, i + 1) - bc (u1, N, i - 1)),
                   2.0 * std::abs (bc (u1, N, i + 1) - bc (u1, N, i)));
    phiL = (bc (u1, N, i) - bc (u1, N, i - 1)) * (bc (u1, N, i - 1) - bc (u1, N, i - 2));
    phiR = (bc (u1, N, i + 1) - bc (u1, N, i)) * (bc (u1, N, i) - bc (u1, N, i - 1));
    VLDeltaL = (phiL > 0) ? (std::copysign (thetaL, (bc (u1, N, i) - bc (u1, N, i - 2)))) : 0;
    VLDeltaR = (phiR > 0) ? (std::copysign (thetaR, (bc (u1, N, i + 1) - bc (u1, N, i + 1)))) : 0;
    fluxL = bc (u1, N, i - 1) + 0.5 * (1.0 - sigma) * VLDeltaL;
    fluxR = bc (u1, N, i) + 0.5 * (1.0 - sigma) * VLDeltaR;
    u[i] = bc(u1, N, i) + sigma * (fluxL - fluxR);
  }
}
