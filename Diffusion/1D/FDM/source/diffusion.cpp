#include <diffusion.h>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Cholesky>
#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;
using namespace std;


void forwardEuler (Ref<VectorXd> u,
                   const int N,
                   const double kappa,
                   const double delta_t) {
  const double h =  1.0 / (N + 1);
  const double mu = delta_t * kappa / (h * h);

  if (mu > 0.5) cerr << "Warning! CFL condition mu = " << mu << " > 0.5; Numerical stability not guaranteed!" << endl;

  static MatrixXd A;
  static MatrixXd B;
  static int oldN;
  if (oldN != N) {
    A = MatrixXd::Zero (N, N);
    A.diagonal (0) = VectorXd::Constant (N, -2);
    A.diagonal (-1) = A.diagonal (1) = VectorXd::Constant (N - 1, 1);
    B = MatrixXd::Identity (N, N) + (mu * A);
    oldN = N;
  }

  u = B * u;
}

void backwardEuler (Ref<VectorXd> u,
                    const int N,
                    const double kappa,
                    const double delta_t) {
  const double h  = 1.0 / (N + 1);
  const double mu = delta_t * kappa / (h * h);

  static MatrixXd A;
  static MatrixXd B;
  static BiCGSTAB<MatrixXd> solver;
  static int oldN;
  if (oldN != N) {
    A = MatrixXd::Zero (N, N);
    A.diagonal (0) = VectorXd::Constant (N, -2);
    A.diagonal (-1) = A.diagonal (1) = VectorXd::Constant (N - 1, 1);
    B = MatrixXd::Identity (N, N) - (mu * A);
    solver.compute (B);
    oldN = N;
  }
  
  VectorXd rhs = u;
  u = solver.solveWithGuess (rhs, u);
}

void crankNicolson (Ref<VectorXd> u,
                    const int N,
                    const double kappa,
                    const double delta_t,
                    const double c) {
  const double h  = 1.0 / (N + 1);
  const double mu = delta_t * kappa / (h * h);

  if (mu > 0.5) cerr << "==> Warning! CFL condition mu = " << mu << " > 0.5; Numerical accuracy not guaranteed!" << endl;

  static MatrixXd A;
  static int oldN;
  if (oldN != N) {
    A = MatrixXd::Zero (N, N);
    A.diagonal (0) = VectorXd::Constant (N, -2);
    A.diagonal (-1) = A.diagonal (1) = VectorXd::Constant (N - 1, 1);
    oldN = N;
  }

  VectorXd rhs = (MatrixXd::Identity (N, N) + c * mu * A) * u;
  MatrixXd B   = (MatrixXd::Identity (N, N) + (c - 1) * mu * A);

  u = B.llt().solve(rhs);
}

void BDF2 (Ref<VectorXd> u,
           Ref<VectorXd> u1,
           const int N,
           const double kappa,
           const double delta_t) {
  const double h  = 1.0 / (N + 1);
  const double mu = delta_t * kappa / (h * h);

  if (mu > 0.5) cerr << "==> Warning! CFL condition mu = " << mu << " > 0.5; Numerical accuracy not guaranteed!" << endl;

  MatrixXd A;
  static MatrixXd B;
  static int oldN;
  if (oldN != N) {
    A = MatrixXd::Zero (N, N);
    A.diagonal (0) = VectorXd::Constant (N, -2);
    A.diagonal (1) = A.diagonal (-1) = VectorXd::Constant (N - 1, 1);
    B = (MatrixXd::Identity(N, N) - 2.0 / 3.0 * mu * A);
    oldN = N;
  }

  VectorXd rhs (N); rhs = u; rhs *= 4.0; rhs /= 3.0;
  VectorXd rhs1 (N); rhs1 = u1; rhs1 /= 3.0; rhs -= rhs1;

  u = B.llt().solve(rhs);

}
