#include <initial-conditions.h>
#include <diffusion.h>
#include <utility.h>
#include <output.h>

#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;
using namespace std;


int main() {
  const int N          = 64;
  const double delta_t = 5e-2;
  const double kappa   = 0.0675;

  VectorXd u (N * N);
  squareWave (u, N);

  backwardEuler (u, N, kappa, delta_t);

  displayField (u, N);

  return 0;
}
