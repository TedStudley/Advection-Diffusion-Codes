#include <initial-conditions.h>
#include <diffusion.h>
#include <output.h>

#include <Eigen/Dense>

#include <fstream>
#include <iostream>

using namespace Eigen;
using namespace std;


int main() {
  const double t0      = 2.5e-5;
  const double mu      = 1;
  const double kappa   = 1.0;
  double t;

  for (int N = 8; N <= 256; N*=2) {
    VectorXd u (N);
    VectorXd u1 (N);
    VectorXd utemp (N);

    double h  = 1.0 / (N + 1);
    double delta_t = mu * h * h / kappa;

    squareWave (u1, N);
    fourierSquare (u, N, kappa, t0);

    t = t0;

    for (int i = 0; t < 0.25; ++i) {
      u = utemp;
      BDF2 (u, u1, N, kappa, delta_t);
      u1 = utemp;
      t += delta_t;
    }

    fourierSquare (u1, N, kappa, t);
    
    VectorXd error = (u - u1);

    cout << N << " " << error.cwiseAbs().maxCoeff() 
              << " " << (error.norm() / error.cols()) << endl;

  }

  return 0;
}
