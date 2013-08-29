#include <initial-conditions.h>
#include <diffusion.h>

#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;
using namespace std;


int main() {
  const int N          = 1024;
  const double delta_t = 5e-6;
  const double kappa   = 0.0675;

  VectorXd u (N);
  squareWave (u, N);
  for (int i = 0; i < 100 ; ++i)
    crankNicolson (u, N, kappa, delta_t);

  cout << u.transpose() << endl;
  
  return 0;
}
