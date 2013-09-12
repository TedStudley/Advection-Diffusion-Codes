#include <cmath>

#include <Eigen/Dense>

using namespace Eigen;


void squareWave (Ref<VectorXd> u) {
  const int N = u.rows();
  const double h = 1.0 / (N + 1);
  double x = h;

  for (int i = 0; i < N; ++i) {
    u[i] = (0.25 < x && x <= 0.75) ? 1 : 0;
    x += h;
  }
}

void fourierSquare (Ref<VectorXd> u,
                    const double kappa,
                    const double t0) {
  const int N = u.rows();
  const double h = 1.0 / (N + 1);
  double x = h;

  int M = 300 / kappa;
  if (N > M) M = N;

  VectorXd bk (M);  
  for (int k = 0; k < M; ++k)
    bk[k] = 2.0 * (cos ((k + 1) * M_PI * 0.25) - cos ((k + 1) * M_PI * 0.75)) / ((k + 1) * M_PI) * exp (-(k + 1) * (k + 1) * kappa * t0 * M_PI * M_PI);
  for (int i = 0; i < N; ++i){
    u[i] = 0;
    for (int k = 0; k < M; ++k) 
      u[i] += bk[k] * sin ((k + 1) * M_PI * x);
    x += h;
  }
}

void sineWave (Ref<VectorXd> u,
               const int k) {
  const int N = u.rows();
  const double h = 1.0 / (N + 1);
  double x = h;

  for (int i = 0; i < N; ++i) {
    u[i] = sin (k * M_PI * x);
    x += h;
  }
}

void sineWave (Ref<VectorXd> u,
               const int k,
               const double kappa,
               const double t0) {
  const int N = u.rows();
  const double h = 1.0 / (N + 1);
  double x = h;

  for (int i = 0; i < N; ++i) {
    u[i] = std::exp(-M_PI * M_PI * k * k * kappa * t0) * std::sin (k * M_PI * x);
    x += h;
  }
}
