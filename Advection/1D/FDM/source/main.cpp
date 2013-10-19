#include <initial-conditions.h>
#include <advection.h>
#include <output.h>

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

int main() {
  const int N          = 256;
  const double delta_t = 0.5 / N;
  const double v       = 1.0;
  double t             = 0; 

  VectorXd u (N);
  VectorXd u1 (N);
  VectorXd uExact (N);
  squareWave (u, N);
  squareWave (u1, N);
  squareWave (uExact, N);

  for (int timestep = 0; t < 1.0; ++timestep) {
    frommVanLeer (u, N, v, delta_t);
    frommMethod (u1, N, v, delta_t);
    t += delta_t;
  }

  displayField (u, N);
  displayField (u1, N);
  displayField (uExact, N);

  return 0;
}
