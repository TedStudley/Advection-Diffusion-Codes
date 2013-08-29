#include <advection.h>
#include <initial-conditions.h>

#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;
using namespace std;

int main() {
  const unsigned int N = 1024;
  const double t0      = 0.00001;
  const double delta_t = 0.0001;
  const double v       = 1.0;

  VectorXd u (N);

  fourierSquare (u, N, t0);

  frommVanLeer (u, N, v, delta_t);

  cout << u.transpose() << endl;

  return 0;
}
