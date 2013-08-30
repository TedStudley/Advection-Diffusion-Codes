#include <cmath>

#include <Eigen/Dense>

using namespace Eigen;


void squareWave (Ref<VectorXd> u,
                 const int N) {
  const double h = 1.0 / (N + 1);
  double x = 0.5 * h;

  for (int i = 0; i < int (N); ++i) {
    u[i] = (0.25 < x && x <= 0.75) ? 1 : 0;
    x += h;
  }
}

void fourierSquare (Ref<VectorXd> u,
                    const int N,
                    const double t0) {
  const double h = 1.0 / (N + 1);
  double x = 0.5 * h;

  unsigned int M = N / 2;

  VectorXd bk (M);  
  for (int k = 0; k < int (M); ++k)
    bk[k] = 2.0 * (cos ((k + 1) * M_PI * 0.25) - cos ((k + 1) * M_PI * 0.75)) / ((k + 1) * M_PI);
  for (int i = 0; i < int (N); ++i){
    u[i] = 0;
    for (int k = 0; k < int (M); ++k) 
      u[i] += bk[k] * exp (-(k + 1) * (k + 1) * t0 * M_PI) * sin ((k + 1) * M_PI * x);
    x += h;
  }
}

