#include <advection.h>
#include <utility.h>

#include <iostream>
#include <cmath>

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void upwindMethod (Ref<VectorXd> u,
                   const int N,
                   const Vector2d v,
                   const double delta_t) {
  const double h     = 1.0 / (N + 1),
               sigma_x = v[0] * delta_t / h,
               sigma_y = v[1] * delta_t / h,
               sigma   = max (sigma_x, sigma_y);

  if (sigma > 1) 
    cerr << "Warning! CFL condition " << sigma << " > 1. Numerical accuracy is not guaranteed." << endl;
  
  VectorXd u1 = u;
  for (int i = 0; i < N; ++i) 
    for (int j = 0; j < N; ++j)
      u[i * N + j] = bc(u1, N, i, j) + sigma_x * (bc(u1, N, i - 1, j) - bc (u1, N, i, j)) + sigma_y * (bc (u1, N, i, j - 1) - bc (u1, N, i, j));
}

void frommMethod (Ref<VectorXd> u,
                  const int N,
                  const Vector2d v,
                  const double delta_t) {
  const double h       = 1.0 / (N + 1),
               sigma_x = v[0] * delta_t / h,
               sigma_y = v[1] * delta_t / h,
               sigma   = max (sigma_x, sigma_y);

  if (sigma > 1)
    cerr << "Warning! CFL condition " << sigma << " > 1. Numerical stability is not guaranteed." << endl;

  VectorXd u1 = u;
  for (int i = 0; i < N; ++i) 
    for (int j = 0; j < N; ++j) 
      u[i * N + j] = bc (u1, N, i, j) + sigma_x * 
                           ((bc (u1, N, i - 1, j) + (1.0 - sigma_x) * 0.25 * (bc (u1, N, i, j) - bc (u1, N, i - 2, j))) - 
                            (bc (u1, N, i, j) + (1.0 - sigma_x) * 0.25 * (bc (u1, N, i + 1, j) - bc (u1, N, i - 1, j)))) +
                                        sigma_y *
                           ((bc (u1, N, i, j - 1) + (1.0 - sigma_y) * 0.25 * (bc (u1, N, i, j) - bc (u1, N, i, j - 2))) -
                            (bc (u1, N, i, j) + (1.0 - sigma_y) * 0.25 * (bc (u1, N, i, j + 1) - bc (u1, N, i, j - 1))));
}

void frommVanLeer (Ref<VectorXd> u,
                   const int N,
                   const Vector2d v,
                   const double delta_t) {
  const double h       = 1.0 / (N + 1),
               sigma_x = v[0] * delta_t / h,
               sigma_y = v[1] * delta_t / h,
               sigma   = max (sigma_x, sigma_y);

  if (sigma > 1)
    cerr << "Warning! CFL conditin " << sigma << " > 1. Numerical stability is not guaranted." << endl;

  VectorXd u1 = u;
