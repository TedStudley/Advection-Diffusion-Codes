#include <diffusion.h>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Cholesky>
#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;
using namespace std;


void forwardEuler (Ref<VectorXd> u,
                   const double delta_t,
                   const double h,
                   const double kappa) {
  const int N = u.rows();
  static MatrixXd B;
  static int oldN;
  if (oldN != N) {
    double mu = delta_t * kappa / (h * h);
    MatrixXd A = MatrixXd::Zero (N, N);
    A.diagonal (0) = VectorXd::Constant (N, 2);
    A.diagonal (-1) = A.diagonal (1) = VectorXd::Constant (N - 1, -1);
    B = MatrixXd::Identity (N, N) - (mu * A);
    oldN = N;
  }

  u = B * u;
}

void backwardEuler (Ref<VectorXd> u,
                    const double delta_t,
                    const double h,
                    const double kappa) {
  const int N = u.rows();
  static MatrixXd B;
  static LDLT<MatrixXd> Bdecomp;
  static int oldN;
  if (oldN != N) {
    double mu = delta_t * kappa / (h * h);
    MatrixXd A = MatrixXd::Zero (N, N);
    A.diagonal (0) = VectorXd::Constant (N, 2);
    A.diagonal (-1) = A.diagonal (1) = VectorXd::Constant (N - 1, -1);
    B = MatrixXd::Identity (N, N) + (mu * A);
    Bdecomp.compute(B);
    oldN = N;
  }
  
  VectorXd rhs = u;
  u = Bdecomp.solve(rhs);
}

void crankNicolson (Ref<VectorXd> u,
                    const double delta_t,
                    const double h,
                    const double kappa) {
  const int N = u.rows();
  
  static MatrixXd B;
  static MatrixXd A;
  static LDLT<MatrixXd> Bdecomp;
  static double mu;
  static int oldN;
  if (oldN != N) {
    mu = delta_t * kappa / (h * h);
    A = MatrixXd::Zero (N, N);
    A.diagonal (0) = VectorXd::Constant (N, 2);
    A.diagonal (-1) = A.diagonal (1) = VectorXd::Constant (N - 1, -1);
    B = (MatrixXd::Identity (N, N) + 0.5 * mu * A);
    Bdecomp.compute(B);
    oldN = N;
  }

  VectorXd rhs = (MatrixXd::Identity (N, N) - 0.5 * mu * A) * u;
  
  u = Bdecomp.solve(rhs);
}

void BDF2 (Ref<VectorXd> u,
           Ref<VectorXd> u1,
           const double delta_t,
           const double h,
           const double kappa) {
  const int N = u.rows();

  static MatrixXd B;
  static LDLT<MatrixXd> Bdecomp;
  static int oldN;
  if (oldN != N) {
    double mu = delta_t * kappa / (h * h);
    MatrixXd A = MatrixXd::Zero (N, N);
    A.diagonal (0) = VectorXd::Constant (N, 2);
    A.diagonal (1) = A.diagonal (-1) = VectorXd::Constant (N - 1, -1);
    B = (MatrixXd::Identity(N, N) + 2.0 / 3.0 * mu * A);
    Bdecomp.compute(B);
    oldN = N;
  }

  VectorXd rhs (N); rhs = u; rhs *= 4.0; rhs /= 3.0;
  VectorXd rhs1 (N); rhs1 = u1; rhs1 /= 3.0; rhs -= rhs1;

  u = Bdecomp.solve (rhs);
}

void TR_BDF2 (Ref<VectorXd> u,
              const double delta_t,
              const double h,
              const double kappa) {
  const int N = u.rows();

  static MatrixXd A;
  static LDLT<MatrixXd> Bsolve;
  static LDLT<MatrixXd> Csolve;
  static double mu;
  static int oldN;
  if (oldN != N) {
    mu = delta_t * kappa / (h * h);
    A = MatrixXd::Zero (N, N);
    A.diagonal (0) = VectorXd::Constant (N, 2);
    A.diagonal (1) = A.diagonal (-1) = VectorXd::Constant (N - 1, -1);
    MatrixXd B = (MatrixXd::Identity(N, N) + delta_t / 4.0 * mu * A);
    MatrixXd C = (MatrixXd::Identity(N, N) + delta_t / 3.0 * mu * A);
    Bsolve.compute(B);
    Csolve.compute(C);
    oldN = N;
  }


  VectorXd rhs (N);

  rhs = (MatrixXd::Identity(N, N) - (delta_t / 4.0 * mu * A)) * u;

  VectorXd u1 (N);
  u1 = Bsolve.solve(rhs);

  rhs = 4.0 / 3.0 * u1 - 1.0 / 3.0 * u;
  u = Csolve.solve(rhs);
}

  

