#include <initial-conditions.h>
#include <diffusion.h>
#include <utility.h>
#include <output.h>

#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;
using namespace std;


int main() {
  const double mu    = 0.30;
  const double kappa = 1.0;
  const int N        = 64;
  const double h     = 1.0 / (N + 1);
  double delta_t = mu * h / kappa;

  VectorXd u (N * N);
  VectorXd u1 (N * N);
  squareWave (u1);
  fourierSquare (u, kappa, delta_t);



  crankNicolson (u, delta_t, h, kappa);
  displayField (u);

  return 0;
}
