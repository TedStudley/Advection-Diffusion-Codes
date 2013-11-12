#pragma once

#define max(a, b) ((a > b) ? a : b)

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <vector>

using namespace Eigen;
using namespace std;

// Forward Grad x (forward-differenced grad)
SparseMatrix<double> fGradX (const int N, const int off = 0) {
  vector<Triplet<double>> tripletList;
  tripletList.reserve (2 * N * N);
  int jb, jf;

  for (int i = 0; i < N; ++i) {
    jb = off; jf = 1 + off;
    while (jb < 0 || jb >= N) (jb = (jb + N) % N);
    while (jf < 0 || jf >= N) (jf = (jf + N) % N);  
    for (int j = 0; j < N; ++j) {
      tripletList.push_back (Triplet<double> (i * N + j, i * N + jb++, -1.0));
      tripletList.push_back (Triplet<double> (i * N + j, i * N + jf++,  1.0));
      if (jb >= N) jb -= N; if (jf >= N) jf -= N;
    }
  }

  SparseMatrix<double> fgradX (N * N, N * N);
  fgradX.setFromTriplets (tripletList.begin (), tripletList.end ());
  return fgradX;
}

// Centered Grad x (center-differenced grad)
SparseMatrix<double> cGradX (const int N, const int off = 0) {
  vector<Triplet<double>> tripletList;
  tripletList.reserve (2 * N * N);
  int jb, jf;

  for (int i = 0; i < N; ++i) {
    jb = N - 1 + off; jf = 1 + off;
    while (jb < 0 || jb >= N) (jb = (jb + N) % N);
    while (jf < 0 || jf >= N) (jf = (jf + N) % N);
    for (int j = 0; j < N; ++j) {
      tripletList.push_back (Triplet<double> (i * N + j, i * N + jb++, -1.0));
      tripletList.push_back (Triplet<double> (i * N + j, i * N + jf++,  1.0));
      if (jb >= N) jb -= N; if (jf >= N) jf -= N;
    }
  }
  
  SparseMatrix<double> cgradX (N * N, N * N);
  cgradX.setFromTriplets (tripletList.begin (), tripletList.end ());

  return cgradX;
}

// Backward Grad x (backward-differenced grad)
SparseMatrix<double> bGradX (const int N, const int off = 0) {
  vector<Triplet<double>> tripletList;
  tripletList.reserve (2 * N * N);
  int jb, jf;

  for (int i = 0; i < N; ++i) {
    jb = N - 1 + off; jf = off;
    while (jb < 0 || jb >= N) (jb = (jb + N) % N);
    while (jf < 0 || jf >= N) (jf = (jf + N) % N);
    for (int j = 0; j < N; ++j) {
      tripletList.push_back (Triplet<double> (i * N + j, i * N + jb++, -1.0));
      tripletList.push_back (Triplet<double> (i * N + j, i * N + jf++,  1.0));
      if (jb >= N) jb -= N; if (jf >= N) jf -= N;
    }
  }
  
  SparseMatrix<double> bgradX (N * N, N * N);
  bgradX.setFromTriplets (tripletList.begin (), tripletList.end ());

  return bgradX;
}

// Forward Grad y (forward-differenced grad)
SparseMatrix<double> fGradY (const int N, const int off = 0) {
  vector<Triplet<double>> tripletList;
  tripletList.reserve (2 * N * N);
  int jb, jf;

  jb = off; jf = N + off;
  while (jb < 0 || jb >= N * N) (jb = (jb + N * N) % (N * N));
  while (jf < 0 || jf >= N * N) (jf = (jf + N * N) % (N * N));
  for (int i = 0; i < N * N; ++i) {
    tripletList.push_back (Triplet<double> (i, jb++, -1.0));
    tripletList.push_back (Triplet<double> (i, jf++,  1.0));
    if (jb >= N * N) jb -= N * N; if (jf >= N * N) jf -= N * N;
  }

  SparseMatrix<double> fgradY (N * N, N * N);
  fgradY.setFromTriplets (tripletList.begin (), tripletList.end ());

  return fgradY;
}

// Centered Grad y (centered-difference grad)
SparseMatrix<double> cGradY (const int N, const int off = 0) {
  vector<Triplet<double>> tripletList;
  tripletList.reserve (2 * N * N);
  int jb, jf;

  jb = N * N - N + off; jf = N + off;
  while (jb < 0 || jb >= N * N) (jb = (jb + N * N) % (N * N));
  while (jf < 0 || jf >= N * N) (jf = (jf + N * N) % (N * N));
  for (int i = 0; i < N * N; ++i) {
    tripletList.push_back (Triplet<double> (i, jb++, -1.0));
    tripletList.push_back (Triplet<double> (i, jf++,  1.0));
    if (jb >= N * N) jb -= N * N; if (jf >= N * N) jf -= N * N;
  }

  SparseMatrix<double> cgradY (N * N, N * N);
  cgradY.setFromTriplets (tripletList.begin (), tripletList.end ());
  
  return cgradY;
}

// Backward Grad y (backward-difference grad)
SparseMatrix<double> bGradY (const int N, const int off = 0) {
  vector <Triplet<double>> tripletList;
  tripletList.reserve (2 * N * N);
  int jb, jf;

  jb = N * N - N + off; jf = off;
  while (jb < 0 || jb >= N * N) (jb = (jb + N * N) % (N * N));
  while (jf < 0 || jf >= N * N) (jf = (jf + N * N) % (N * N));
  for (int i = 0; i < N * N; ++i) {
    tripletList.push_back (Triplet<double> (i, jb++, -1.0));
    tripletList.push_back (Triplet<double> (i, jf++,  1.0));
    if (jb >= N * N) jb -= N * N; if (jf >= N * N) jf -= N * N;
  }

  SparseMatrix<double> bgradY (N * N, N * N);
  bgradY.setFromTriplets (tripletList.begin (), tripletList.end ());

  return bgradY;
}
