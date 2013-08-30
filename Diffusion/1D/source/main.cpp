#include <initial-conditions.h>
#include <diffusion.h>
#include <output.h>

#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;
using namespace std;


int main() {
  const int N          = 1024;
  const double t0      = 1e-7;
  const double delta_t = 5e-6;
  const double kappa   = 0.0675;

  VectorXd u (N);
  fourierSquare (u, N, t0);
  for (int i = 0; i < 1 ; ++i)
    crankNicolson (u, N, kappa, delta_t);

  displayField (u, N);

  return 0;
}
