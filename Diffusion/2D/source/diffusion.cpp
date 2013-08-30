#include <diffusion.h>

#include <Eigen/IterativeLinearSolvers>
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
    MatrixXd BBlock = MatrixXd::Zero (N, N);
    BBlock.diagonal (0) = VectorXd::Constant (N, -4);
    BBlock.diagonal (1) = BBlock.diagonal (-1) = VectorXd::Constant (N - 1, 1);
    A = MatrixXd::Zero (N * N,  N * N);
    A.block(0, 0, N, N) = BBlock;
    for (int i = 0; i < N - 1; ++i) {
      A.block(N * (i + 1), N * (i + 1), N, N) = BBlock;
      A.block(N * (i + 1), N * i, N, N) = A.block(N * i, N * (i + 1), N, N) =
        MatrixXd::Identity(N, N);
    }
    B = MatrixXd::Identity (N * N, N * N) + (mu * A);
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
  static int oldN;
  static BiCGSTAB<MatrixXd> solver;
  if (oldN != N) {
    MatrixXd BBlock = MatrixXd::Zero (N, N);
    BBlock.diagonal (0) = VectorXd::Constant (N, -4);
    BBlock.diagonal (1) = BBlock.diagonal (-1) = VectorXd::Constant (N - 1, 1);
    A = MatrixXd::Zero (N * N, N * N);
    A.block(0, 0, N, N) = BBlock;
    for (int i = 0; i < N - 1; ++i) {
      A.block(N * (i + 1), N * (i + 1), N, N) = BBlock;
      A.block(N * (i + 1), N * i, N, N) = A.block(N * i, N * (i + 1), N, N) =
        MatrixXd::Identity(N, N);
    }
    
    B = MatrixXd::Identity (N * N, N * N) - (mu * A);
    solver.compute (B);
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
  static BiCGSTAB<MatrixXd> solver;
  static int oldN;
  if (oldN != N) {
    MatrixXd BBlock = MatrixXd::Zero (N, N);
    BBlock.diagonal (0) = VectorXd::Constant (N, -4);
    BBlock.diagonal (1) = BBlock.diagonal (-1) = VectorXd::Constant (N - 1, 1);
    A = MatrixXd::Zero (N * N, N * N);
    A.block(0, 0, N, N) = BBlock;
    for (int i = 0; i < N - 1; ++i) {
      A.block (N * (i + 1), N * (i + 1), N, N) = BBlock;
      A.block (N * (i + 1), N * i, N, N) = A.block (N * i, N * (i + 1), N , N) =
        MatrixXd::Identity (N, N);
    }
  }

  VectorXd rhs = (MatrixXd::Identity (N * N, N * N) + c * mu * A) * u;
  MatrixXd B   = (MatrixXd::Identity (N * N, N * N) + (c - 1) * mu * A);

  solver.compute(B);

  u = solver.solveWithGuess(rhs, u);
}
