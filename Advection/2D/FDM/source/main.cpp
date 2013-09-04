#include <initial-conditions.h>
#include <advection.h>
#include <output.h>

#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;
using namespace std;

int main() {
  const int N          = 64;
  const double delta_t = 0.001;
  Vector2d v           = {1.0, 1.0};

  VectorXd u (N * N);

  squareWave (u, N);

  for (int i = 0; i < 1; ++i) {
    cerr << "timestep " << i << ":" << endl;
    frommVanLeer (u, N, v, delta_t);
  }
  
  displayField (u, N);

  return 0;
}
