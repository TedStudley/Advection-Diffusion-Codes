#include <cmath>

#include <Eigen/Dense>

using namespace Eigen;


void squareWave (Ref<VectorXd> u) {
  const int N = sqrt (u.rows ());
  const double h = 1.0 / (N + 1);

  Vector2d x = Vector2d::Constant (0.5 * h);

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      u[i * N + j] = 
        ((0.25 < x(0) && x(0) <= 0.75) &&
         (0.25 < x(1) && x(1) <= 0.75)) ? 1 : 0;
      x(0) += h;
    }
    x(0) = 0.5 * h; x(1) += h;
  }
}

void fourierSquare (Ref<VectorXd> u,
                    const double kappa,
                    const double t0) {
  const int N = sqrt (u.rows ());
  const double h = 1.0 / (N + 1);
  Vector2d x = Vector2d::Constant (0.5 * h);

  VectorXd bk (N);  
  for (int k = 0; k < N; ++k)
    bk[k] = 2.0 * (cos ((k + 1) * M_PI * 0.25) - cos ((k + 1) * M_PI * 0.75)) / ((k + 1) * M_PI);
  for (int i = 0; i < N; ++i){
    for (int j = 0; j < N; ++j) {
      u[i * N + j] = 0;
      for (int k1 = 0; k1 < N; ++k1) {
        for (int k2 = 0; k2 < N; ++k2) {
          u[i * N + j] += bk[k1] * bk[k2] * exp (-(k1 + 1) * (k1 + 1) * kappa * t0 * M_PI * M_PI) * exp (-(k2 + 1) * (k2 + 1) * t0 * M_PI * M_PI) * sin ((k1 + 1) * M_PI * x[0]) * sin ((k2 + 1) * M_PI * x[1]);
        }
      }
      x[0] += h;
    }
    x[0] = 0.5 * h; x[1] += h;
  }
}
