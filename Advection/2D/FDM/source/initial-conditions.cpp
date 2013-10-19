#include <cmath>

#include <Eigen/Dense>

using namespace Eigen;


void squareWave (Ref<VectorXd> u) {
  const int N = sqrt (u.rows());
  const double h = 1.0 / (N + 1);
  Vector2d x = Vector2d::Constant (h);

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      u[i * N + j] = 
        ((0.25 < x(0) && x(0) <= 0.75) &&
         (0.25 < x(1) && x(1) <= 0.75)) ? 1 : 0;
      x(0) += h;
    }
    x(0) = h; x(1) += h;
  }
}

void fourierSquare (Ref<VectorXd> u,
                    const double kappa,
                    const double t0) {
  const int N = sqrt (u.rows ());
  const double h = 1.0 / (N + 1);
  Vector2d x = Vector2d::Constant (h);

  int M = 300;
  if (N < M) M = N;

  VectorXd bk (M);  
  for (int k = 0; k < M; ++k)
    bk[k] = 2.0 * (cos ((k + 1) * M_PI * 0.25) - cos ((k + 1) * M_PI * 0.75)) / ((k + 1) * M_PI) * exp (-(k + 1) * (k + 1) * kappa * t0 * M_PI * M_PI);
  for (int i = 0; i < N; ++i){
    for (int j = 0; j < N; ++j) {
      u[i * N + j] = 0;
      for (int k1 = 0; k1 < M; ++k1) {
        for (int k2 = 0; k2 < M; ++k2) {
          u[i * N + j] += bk[k1] * bk[k2] * sin ((k1 + 1) * M_PI * x[0]) * sin ((k2 + 1) * M_PI * x[1]);
        }
      }
      x[0] += h;
    }
    x[0] = h; x[1] += h;
  }
}

void sineWave (Ref<VectorXd> u,
               const int k1,
               const int k2) {
  const int N = sqrt (u.rows ());
  const double h = 1.0 / (N + 1);
  Vector2d x = Vector2d::Constant (h);

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      u[i * N + j] = sin (k1 * M_PI * x[0]) * sin (k1 * M_PI * x[1]);
      x[0] += h;
    }
    x[0] = h; x[1] += h;
  }
}

void sineWave (Ref<VectorXd> u,
               const int k1, 
               const int k2,
               const double kappa,
               const double t0) {
  const int N = sqrt (u.rows ());
  const double h = 1.0 / (N + 1);
  Vector2d x = Vector2d::Constant (h);

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      u[i * N + j] = sin (k1 * M_PI * x[0]) * sin (k1 * M_PI * x[1]) * exp(-k1 * k1 * kappa * t0 * M_PI * M_PI) * exp(-k2 * k2 * kappa * t0 * M_PI * M_PI);
      x[0] += h;
    }
    x[0] = h; x[1] += h;
  }
}
