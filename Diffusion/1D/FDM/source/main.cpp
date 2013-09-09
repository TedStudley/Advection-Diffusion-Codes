#include <initial-conditions.h>
#include <diffusion.h>
#include <output.h>
#include <norms.h>

#include <Eigen/Dense>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace Eigen;
using namespace std;


int main() {
  const double mu      = 0.5;
  const double kappa   = 1.0;
  const double t0      = 0.0001;
  const double T       = 0.1;
  const int k = 1;
  double t;

  for (int N = 16; N <= 8192; N*=2) {
    VectorXd u (N);
    VectorXd u1 (N);
    VectorXd utemp (N);

    double h  = 1.0 / (N + 1);
    double delta_t = mu * h / kappa;

    int N_timestep = T / delta_t;
    if (abs (N_timestep * delta_t - T) > 10e-16) {
      delta_t = T / (++N_timestep);
    }

    fourierSquare (u1, kappa, delta_t);
    fourierSquare (u, kappa, 2 * delta_t);

    t = 2 * delta_t;

    for (int i = 0; t < T; ++i) {
      utemp  = u;
      BDF2 (u, u1, delta_t, h, kappa);
      u1 = utemp;
      t += delta_t;
    }

    fourierSquare (u1, kappa, t);

    VectorXd error = (u - u1);

    cerr << N << " " << maxNorm (u - u1) 
         << " " << sqrt(error.array().square().sum()) / N << endl;

  }

  return 0;
}
