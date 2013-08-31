#include <initial-conditions.h>
#include <diffusion.h>
#include <output.h>

#include <Eigen/Dense>

#include <fstream>
#include <iostream>

using namespace Eigen;
using namespace std;


int main() {
  const int N          = 1024;
  const double t0      = 2.5e-5;
  const double delta_t = 5e-4;
  const double kappa   = 1.0;
  VectorXd utemp (N);
  double t;

  VectorXd u (N);
  VectorXd u1 (N);
  fourierSquare (u1, N, kappa, t0);
  fourierSquare (u, N, kappa, delta_t);

  ofstream initialOut ("initial.dat");
  displayField (u, N, initialOut);

  t = delta_t;

  for (int i = 0; i < 1; ++i) {
    utemp = u;
    crankNicolson (u, N, kappa, delta_t);
    //BDF2 (u, u1, N, kappa, delta_t);
    u1 = utemp;
    t += delta_t;
  }

  fourierSquare (u1, N, kappa, t);

  ofstream calculatedOut ("log");
  ofstream exactOut ("exact.dat");

  displayField (u, N, calculatedOut);
  displayField (u1, N, exactOut);

  return 0;
}
