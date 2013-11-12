#include <advection.h>
#include <matforms.h>
#include <utility.h>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <Eigen/Cholesky>

#include <iostream>
#include <vector>
#include <cmath>

using namespace Eigen;
using namespace std;

void upwindMethod (Ref<VectorXd> u,
                   const double dt,
                   const double h,
                   const Vector2d v) {
  const int N = sqrt (u.rows ());
  static int oldN;

  static SparseMatrix<double> gradX, gradY;

  if (oldN != N) {
    gradX = bGradX (N);
    gradY = bGradY (N);

    gradX.makeCompressed (); gradY.makeCompressed ();

    oldN = N;
  }
  
  VectorXd u1 = u;
  u = u1 - (v[0] * dt / (2 * h)) * gradX * u1;
  u1 = u - (v[1] * dt / h) * gradY * u;
  u = u1 - (v[0] * dt / (2 * h)) * gradX * u1;

}

void frommMethod (Ref<VectorXd> u,
                  const double dt,
                  const double h,
                  const Vector2d v) {
  const int N = sqrt (u.rows ());
  static int oldN;

  static SparseMatrix<double> gradX, gradY, sigmaX, sigmaY, temp;

  if (oldN != N) {
    gradX = bGradX (N);
    gradY = bGradY (N);

    sigmaX = cGradX (N) - cGradX (N, -1); 
    sigmaY = cGradY (N) - cGradY (N, -1);

    gradX.makeCompressed (); gradY.makeCompressed ();
    sigmaX.makeCompressed (); sigmaY.makeCompressed ();

    oldN = N;
  }

  VectorXd u1 = u;
  u = u1 - ((v[0] * dt / (2 * h)) * gradX + (v[0] * dt / (8.0 * h)) * (1 - v[0] * dt / (2.0 * h)) * sigmaX) * u1;
  u1 = u - ((v[1] * dt / h) * gradY + (v[1] * dt / (4.0 * h)) * (1 - v[1] * dt / h) * sigmaY) * u;
  u = u1 - ((v[0] * dt / (2 * h)) * gradX + (v[0] * dt / (8.0 * h)) * (1 - v[0] * dt / (2.0 * h)) * sigmaX) * u1;
}

void beamWarming (Ref<VectorXd> u,
                  const double dt,
                  const double n,
                  const VectorXd v) {
  const int N = sqrt (u.rows ());
  static int oldN;

  static SparseMatrix<double> gradX, gradY, sigmaX, sigmaY;

  if (oldN != N) {
    gradX = bGradX (N);
    gradY = bGradY (N);

    sigmaX = bGradX (N) - bGradX (N, -1);
    sigmaY = bgradY (N) - bGradY (N, -1);

    gradX.makeCompressed (); gradY.makeCompressed ();
    sigmaX.makeCompressed (); sigmaY.makeCompressed ();

    oldN = N;
  }

  VectorXd u1 = u;
  
  u = u1 - ((v[0] * dt / (2 * h)) * gradX + (v[0] * dt / (4.0 * h)) * (1 - v[0] * dt / (2.0 * h)) * sigmaX) * u1;
  u1 = u - ((v[1] * dt / h) * gradY + (v[1] * dt / (2.0 * h)) * (1 - v[1] * dt / h) * sigmaY) * u;
  u = u1 - ((v[0] * dt / (2 * h)) * gradX + (v[0] * dt / (4.0 * h)) * (1 - v[0] * dt / (2.0 * h)) * sigmaX) * u1;
}

void laxWendroff (Ref<VectorXd> u,
                  const double dt,
                  const double h,
                  const Vector2d v) {
  const int N = sqrt (u.rows ());
  static int oldN;

  static SparseMatrix<double> gradX, gradY, sigmaX, sigmaY;

  if (oldN != N) {
    gradX = bGradX (N);
    gradY = bGradY (N);

    sigmaX = fGradX (N) - fGradX (N, -1);
    sigmaY = fGradY (N) - fGradY (N, -1);

    gradX.makeCompressed (); gradY.makeCompressed ();
    sigmaX.makeCompressed (); sigmaY.makeCompressed ();
    
    oldN = N;
  }

  VectorXd u1 = 1;

  u = u1 - ((v[0] * dt / (2 * h)) * gradX + (v[0] * dt / (4.0 * h)) * (1 - v[0] * dt / (2.0 * h)) * sigmaX) * u1;
  u1 = u - ((v[1] * dt / h) * gradY + (v[1] * dt / (2.0 * h)) * (1 - v[1] * dt / h) * sigmaY) * u;
  u = u1 - ((v[0] * dt / (2 * h)) * gradX + (v[0] * dt / (4.0 * h)) * (1 - v[0] * dt / (2.0 * h)) * sigmaX) * u1;
}


void frommVanLeer (Ref<VectorXd> u,
                   const double dt,
                   const double h,
                   const Vector2d v) {
  const int N = sqrt (u.rows ());

  double thetaL_x, thetaR_x, phiL_x, phiR_x, VLDeltaL_x, VLDeltaR_x, fluxL_x, fluxR_x;
  double thetaL_y, thetaR_y, phiL_y, phiR_y, VLDeltaL_y, VLDeltaR_y, fluxL_y, fluxR_y;

  VectorXd u1 = u;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      thetaL_x = min3 (2.0 * std::abs (bc (u1, N, i - 1, j) - bc (u1, N, i - 2, j)),
                       0.5 * std::abs (bc (u1, N, i, j) - bc (u1, N, i - 2, j)),
                       2.0 * std::abs (bc (u1, N, i, j) - bc (u1, N, i - 1, j)));
      thetaR_x = min3 (2.0 * std::abs (bc (u1, N, i, j) - bc (u1, N, i - 1, j)),
                       0.5 * std::abs (bc (u1, N, i + 1, j) - bc (u1, N, i - 1, j)),
                       2.0 * std::abs (bc (u1, N, i + 1, j) - bc (u1, N, i, j)));
      
      phiL_x = (bc (u1, N, i, j) - bc (u1, N, i - 1, j)) *
               (bc (u1, N, i - 1, j) - bc (u1, N, i - 2, j));
      phiR_x = (bc (u1, N, i + 1, j) - bc (u1, N, i, j)) *
               (bc (u1, N, i, j) - bc (u1, N, i - 1, j));
      
      VLDeltaL_x = (phiL_x > 0) ? copysign (thetaL_x, (bc (u1, N, i, j) - bc (u1, N, i - 2, j))) : 0;
      VLDeltaR_x = (phiR_x > 0) ? copysign (thetaR_x, (bc (u1, N, i + 1, j) - bc (u1, N, i - 1, j))) : 0;

      fluxL_x = bc (u1, N, i - 1, j) + 0.5 * (1.0 - v[0] * dt / (2.0 * h)) * VLDeltaL_x;
      fluxR_x = bc (u1, N, i, j) + 0.5 * (1.0 - v[0] * dt / (2.0 * h)) * VLDeltaR_x;


      u[i * N + j] = bc (u1, N, i, j) + v[0] * dt / (2.0 * h) * (fluxL_x - fluxR_x);
    }
  }

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      thetaL_y = min3 (2.0 * std::abs (bc (u, N, i, j - 1) - bc (u, N, i, j - 2)),
                       0.5 * std::abs (bc (u, N, i, j) - bc (u, N, i, j - 2)),
                       2.0 * std::abs (bc (u, N, i, j) - bc (u, N, i, j - 1)));
      thetaR_y = min3 (2.0 * std::abs (bc (u, N, i, j) - bc (u, N, i, j - 1)),
                       0.5 * std::abs (bc (u, N, i, j + 1) - bc (u, N, i, j - 1)),
                       2.0 * std::abs (bc (u, N, i, j + 1) - bc (u, N, i, j)));

      phiL_y = (bc (u, N, i, j) - bc (u, N, i, j - 1)) *
               (bc (u, N, i, j - 1) - bc (u, N, i, j - 2));
      phiR_y = (bc (u, N, i, j + 1) - bc (u, N, i, j)) *
               (bc (u, N, i, j) - bc (u, N, i, j - 1));

      VLDeltaL_y = (phiL_y > 0) ?
                    copysign (thetaL_y,
                               (bc (u, N, i, j) - bc (u, N, i, j - 2)))
                     : 0;
      VLDeltaR_y = (phiR_y > 0) ?
                     copysign (thetaR_y,
                               (bc (u, N, i, j + 1) - bc (u, N, i, j - 1)))
                     : 0;

      fluxL_y = bc (u, N, i, j - 1) + 0.5 * (1.0 - v[1] * dt / h) * VLDeltaL_y;
      fluxR_y = bc (u, N, i, j) + 0.5 * (1.0 - v[1] * dt / h) * VLDeltaR_y;

      u1[i * N + j] = bc (u, N, i, j) + v[1] * dt / h * (fluxL_y - fluxR_y);
    }
  }

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      thetaL_x = min3 (2.0 * std::abs (bc (u1, N, i - 1, j) - bc (u1, N, i - 2, j)),
                       0.5 * std::abs (bc (u1, N, i, j) - bc (u1, N, i - 2, j)),
                       2.0 * std::abs (bc (u1, N, i, j) - bc (u1, N, i - 1, j)));
      thetaR_x = min3 (2.0 * std::abs (bc (u1, N, i, j) - bc (u1, N, i - 1, j)),
                       0.5 * std::abs (bc (u1, N, i + 1, j) - bc (u1, N, i - 1, j)),
                       2.0 * std::abs (bc (u1, N, i + 1, j) - bc (u1, N, i, j)));
      
      phiL_x = (bc (u1, N, i, j) - bc (u1, N, i - 1, j)) *
               (bc (u1, N, i - 1, j) - bc (u1, N, i - 2, j));
      phiR_x = (bc (u1, N, i + 1, j) - bc (u1, N, i, j)) *
               (bc (u1, N, i, j) - bc (u1, N, i - 1, j));
      
      VLDeltaL_x = (phiL_x > 0) ? 
                     copysign (thetaL_x, 
                               (bc (u1, N, i, j) - bc (u1, N, i - 2, j)))
                     : 0;
      VLDeltaR_x = (phiR_x > 0) ?
                     copysign (thetaR_x,
                               (bc (u1, N, i + 1, j) - bc (u1, N, i - 1, j)))
                     : 0;

      fluxL_x = bc (u1, N, i - 1, j) + 0.5 * (1.0 - v[0] * dt / (2.0 * h)) * VLDeltaL_x;
      fluxR_x = bc (u1, N, i, j) + 0.5 * (1.0 - v[0] * dt / (2.0 * h)) * VLDeltaR_x;


      u[i * N + j] = bc (u1, N, i, j) + v[0] * dt / (2.0 * h) * (fluxL_x - fluxR_x);
    }
  }
}

