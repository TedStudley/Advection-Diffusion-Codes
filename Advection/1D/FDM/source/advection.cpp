#include <advection.h>
#include <matforms.h>
#include <utility.h>

#include <iostream>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

void upwindMethod (Ref<VectorXd> u,
                   const double dt,
                   const double h,
                   const VectorXd v) {
  const int N = u.rows ();
  static int oldN;

  static SparseMatrix<double> gradX;

  if (oldN != N) {
    gradX = bGradX (N);

    gradX.makeCompressed ();

    oldN = N;
  }

  VectorXd u1 = u;

  u = u1 - (v[0] * dt / h) * gradX * u1;
}

void frommMethod (Ref<VectorXd> u,
                  const double dt,
                  const double h,
                  const VectorXd v) {
  const int N = u.rows ();
  static int oldN;

  static SparseMatrix<double> gradX, sigmaX;

  if (oldN != N) {
    gradX = bGradX (N);
    sigmaX = cGradX (N) - cGradX (N, -1);

    gradX.makeCompressed ();
    sigmaX.makeCompressed ();

    oldN = N;
  }

  VectorXd u1 = u;
  
  u = u1 - ((v[0] * dt / h) * gradX + (v[0] * dt / (4.0 * h) - v[0] * v[0] * dt * dt / (4.0 * h * h)) * sigmaX) * u1;
}

void beamWarming (Ref<VectorXd> u,
                  const double dt,
                  const double h,
                  const VectorXd v) {
  const int N = u.rows ();
  static int oldN;

  static SparseMatrix<double> gradX, sigmaX;

  if (oldN != N) {
    gradX= bGradX (N);
    sigmaX = bGradX (N) - bGradX (N, -1);

    gradX.makeCompressed ();
    sigmaX.makeCompressed ();

    oldN = N;
  }

  VectorXd u1 = u;
  
  u = u1 - ((v[0] * dt / h) * gradX + (v[0] * dt / (2.0 * h) - v[0] * v[0] * dt * dt / (2.0 * h * h)) * sigmaX) * u1;
}

void laxWendroff (Ref<VectorXd> u,
                  const double dt,
                  const double h,
                  const VectorXd v) {
  const int N = u.rows ();
  static int oldN;

  static SparseMatrix<double> gradX, sigmaX;

  if (oldN != N) {
    gradX = bGradX (N);
    sigmaX = fGradX (N) - fGradX (N, -1);
  
    gradX.makeCompressed ();
    sigmaX.makeCompressed ();

    oldN = N;
  }

  VectorXd u1 = u;

  u = u1 - ((v[0] * dt / h) * gradX + (v[0] * dt / (2.0 * h) - v[0] * v[0] * dt * dt / (2.0 * h * h)) * sigmaX) * u1;
}

void frommVanLeer (Ref<VectorXd> u,
                   const double dt,
                   const double h,
                   const VectorXd v) {
  const int N = u.rows ();
/* 
  static int oldN;

  static MatrixXd grad;
  static MatrixXd cDiff, lDiff, rDiff;

  if (oldN != N) {
    grad = MatrixXd::Zero (N, N);
    grad.diagonal () = VectorXd::Constant (N, 1.0);
    grad.diagonal (-1) = VectorXd::Constant (N - 1, -1.0);
    grad (0, N - 1) = -1.0;

    cDiff = MatrixXd::Zero (N, N);
    cDiff.diagonal (1) = VectorXd::Constant (N - 1, 1.0);
    cDiff.diagonal (-1) = VectorXd::Constant (N - 1, -1.0);
    cDiff (0, N - 1) = -1.0; cDiff (N - 1, 0) = 1.0;

    lDiff = MatrixXd::Zero (N, N);
    lDiff.diagonal () = VectorXd::Constant (N, 1.0);
    lDiff.diagonal (-1) = VectorXd::Constant (N - 1, -1.0);
    lDiff (0, N - 1) = -1.0;

    rDiff = MatrixXd::Zero (N, N);
    rDiff.diagonal () = VectorXd::Constant (N, -1.0);
    rDiff.diagonal (1) = VectorXd::Constant (N - 1, 1.0);
    rDiff (N - 1, 0) = 1.0;
  }

  VectorXd theta = (2*lDiff*u).cwiseAbs().cwiseMin((cDiff*u/2.0).cwiseAbs().cwiseMin((2*rDiff*u).cwiseAbs()));

  VectorXd phi = (rDiff*u).cwiseProduct((lDiff*u));

  VectorXd delta = ((cDiff*u).cwiseAbs().cwiseProduct(theta.unaryExpr(unarySign<double>()))).binaryExpr(phi, buildDelta<double>());

  VectorXd u1 = (v[0] * dt / h) * grad * u + (v[0] * dt / (2 * h)) * (1.0 - v[0] * dt / h) * grad * delta;

  u -= u1;
*/
  const double sigma = v[0] * dt / h;

  double thetaL, thetaR, phiL, phiR, VLDeltaL, VLDeltaR, fluxL, fluxR;

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
