#include <advection.h>
#include <initial-conditions.h>

#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;
using namespace std;

int main() {
  const int N          = 1024;
  const double delta_t = 0.001;
  const double v       = 0.5;

  VectorXd u (N);

  squareWave (u, N);

  for (int timestep = 0; timestep < 500; ++timestep)
    frommVanLeer (u, N, v, delta_t);

  cout << u.transpose() << endl;

  return 0;
}
