#include <diffusion.h>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <Eigen/Cholesky>

#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;


void forwardEuler (Ref<VectorXd> u,
                   const double delta_t,
                   const double h,
                   const double kappa) {
  const int N = sqrt (u.rows ());
  static SparseMatrix<double> B;
  static int oldN;
  if (oldN != N) {
    double mu = delta_t * kappa / (h * h);
    if (mu > 0.25)
      cerr << "CAUTION: If the computed solution oscillates, this could be because mu = "
           << mu << " > 0.25" << endl;
    vector<Triplet<double>> tripletList;
    tripletList.reserve (5 * N * N);
    tripletList.push_back (Triplet<double> (0, 0, 1 - 4 * mu));
    for(int i = 1; i < N * N; ++i) {
      tripletList.push_back (Triplet<double> (i, i, 1 - 4 * mu));
      tripletList.push_back (Triplet<double> (i, i - 1, mu));
      tripletList.push_back (Triplet<double> (i - 1, i, mu));
    }
    for (int i = 0; i < N * N - N; ++i) {
      tripletList.push_back (Triplet<double> (i, i + N, mu));
      tripletList.push_back (Triplet<double> (i + N, i, mu));
    }


    B.resize (N * N, N * N);
    B.setFromTriplets (tripletList.begin (), tripletList.end ());

    oldN = N;
  }

  u = B * u;
}

void backwardEuler (Ref<VectorXd> u,
                    const double delta_t,
                    const double h,
                    const double kappa) {
  const int N = sqrt(u.rows());
  static SparseMatrix<double> B;
  static SimplicialLLT<SparseMatrix<double>> Bdecomp;
  static int oldN;
  if (oldN != N) {
    double mu = delta_t * kappa / (h * h);
    vector<Triplet<double>> tripletList;
    tripletList.reserve (5 * N * N);
    tripletList.push_back (Triplet<double>(0, 0, -4 * mu));
    for (int i = 1; i < N * N; ++i) {
      tripletList.push_back (Triplet<double> (i, i, -4 * mu));
      tripletList.push_back (Triplet<double> (i, i - 1, mu));
      tripletList.push_back (Triplet<double> (i - 1, i, mu));
    }
    for (int i = 0; i < N * N - N; ++i) {
      tripletList.push_back (Triplet<double> (i, i + N, mu));
      tripletList.push_back (Triplet<double> (i + N, i, mu));
    }

    SparseMatrix<double> A(N * N, N * N);
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    B.resize (N * N, N * N);
    B.setIdentity();
    B -= A;
    Bdecomp.compute (B);
    oldN = N;
  }
  
  VectorXd rhs = u;
  u = Bdecomp.solve (rhs);
}

void crankNicolson (Ref<VectorXd> u,
                    const double delta_t,
                    const double h,
                    const double kappa) {
  const int N = sqrt(u.rows());
  static SparseMatrix<double> C;
  static SimplicialLLT<SparseMatrix<double>> Bdecomp;
  static double mu;
  static int oldN;
  if (oldN != N) {
    mu = delta_t * kappa / (h * h);
    
    vector<Triplet<double>> tripletList;
    tripletList.reserve (5 * N * N);
    tripletList.push_back (Triplet<double> (0, 0, -4 * mu));
    for (int i = 1; i < N * N; ++i) {
      tripletList.push_back (Triplet<double> (i, i, -4 * mu));
      tripletList.push_back (Triplet<double> (i, i - 1, mu));
      tripletList.push_back (Triplet<double> (i - 1, i, mu));
    }
    for (int i = 0; i < N * N - N; ++i) {
      tripletList.push_back (Triplet<double> (i, i + N, mu));
      tripletList.push_back (Triplet<double> (i + N, i, mu));
    }

    SparseMatrix<double> A (N * N, N * N);
    A.setFromTriplets (tripletList.begin (), tripletList.end ());

    SparseMatrix<double> B (N * N, N * N);
    B.setIdentity();
    B -= 0.5 * A;
    Bdecomp.compute (B);

    C.resize (N * N, N * N);
    C.setIdentity();
    C += 0.5 * A;

    oldN = N;
  }

  VectorXd rhs = C * u;

  u = Bdecomp.solve (rhs);
}

void BDF2 (Ref<VectorXd> u,
           Ref<VectorXd> u1,
           const double delta_t,
           const double h,
           const double kappa) {
  const int N = sqrt(u.rows());

  static SimplicialLDLT<SparseMatrix<double>> BDecomp;

  static int oldN;
  if (oldN != N) {
    double mu = delta_t * kappa / (h * h);
    
    vector<Triplet<double>> tripletList;
    tripletList.reserve (5 * N * N);
    tripletList.push_back (Triplet<double> (0, 0, -4 * mu));
    for (int i = 1; i < N * N; ++i) {
      tripletList.push_back (Triplet<double> (i, i, -4 * mu));
      tripletList.push_back (Triplet<double> (i, i - 1, mu));
      tripletList.push_back (Triplet<double> (i - 1, i, mu));
    }
    for (int i = 0; i < N * N - N; ++i) {
      tripletList.push_back (Triplet<double> (i, i + N, mu));
      tripletList.push_back (Triplet<double> (i + N, i, mu));
    }

    SparseMatrix<double> A (N * N, N * N);
    A.setFromTriplets (tripletList.begin (), tripletList.end ());

    SparseMatrix<double> B (N * N, N * N);
    B.setIdentity ();
    B -= 2.0 / 3.0 * A;
    BDecomp.compute (B);

    oldN = N;
  }

  VectorXd rhs1 (N * N); rhs1 = u; rhs1 *= 4.0; rhs1 /= 3.0;
  VectorXd rhs2 (N * N); rhs2 = u1; rhs2 /= 3.0; 
  VectorXd rhs = rhs1 - rhs2;

  u = BDecomp.solve (rhs);
}
