#pragma once

#include <Eigen/Sparse>

#include <vector>

using namespace Eigen;
using namespace std;

// Forward Grad x (forward-differenced grad)
SparseMatrix<double> fGradX (const int N, const int off = 0) {
  vector<Triplet<double>> tripletList;
  tripletList.reserve (2 * N);
  int jb, jf;

  jb = off; jf = 1 + off;
  while (jb < 0 || jb >= N) (jb = (jb + N) % N);
  while (jf < 0 || jf >= N) (jf = (jf + N) % N);
  for (int i = 0; i < N; ++i) {
    tripletList.push_back (Triplet<double> (i, jb++, -1.0));
    tripletList.push_back (Triplet<double> (i, jf++,  1.0));
    if (jb >= N) jb -= N; if (jf >= N) jf -= N;
  }

  SparseMatrix<double> fgradX (N, N);
  fgradX.setFromTriplets (tripletList.begin (), tripletList.end ());
  return fgradX;
}

// Centered Grad x (center-differenced grad)
SparseMatrix<double> cGradX (const int N, const int off = 0) {
  vector<Triplet<double>> tripletList;
  tripletList.reserve (2 * N);
  int jb, jf;

  jb = N - 1 + off; jf = 1 + off;
  while (jb < 0 || jb >= N) (jb = (jb + N) % N);
  while (jf < 0 || jf >= N) (jf = (jf + N) % N);
  for (int i = 0; i < N; ++i) {
    tripletList.push_back (Triplet<double> (i, jb++, -1.0));
    tripletList.push_back (Triplet<double> (i, jf++,  1.0));
    if (jb >= N) jb -= N; if (jf >= N) jf -= N;
  }

  SparseMatrix<double> cgradX (N, N);
  cgradX.setFromTriplets (tripletList.begin (), tripletList.end ());
  return cgradX;
}

// Backward Grad x (backward-differenced grad)
SparseMatrix<double> bGradX (const int N, const int off = 0) {
  vector<Triplet<double>> tripletList;
  tripletList.reserve (2 * N);
  int jb, jf;

  jb = N - 1 + off; jf = off;
  while (jb < 0 || jb >= N) (jb = (jb + N) % N);
  while (jf < 0 || jf >= N) (jf = (jf + N) % N);
  for (int i = 0; i < N; ++i) {
    tripletList.push_back (Triplet<double> (i, jb++, -1.0));
    tripletList.push_back (Triplet<double> (i, jf++,  1.0));
    if (jb >= N) jb -= N; if (jf >= N) jf -= N;
  }

  SparseMatrix<double> bgradX (N, N);
  bgradX.setFromTriplets (tripletList.begin (), tripletList.end ());
  return bgradX;
}
    
