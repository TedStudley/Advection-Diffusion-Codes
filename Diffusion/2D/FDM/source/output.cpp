#include <output.h>

#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;
using namespace std;

void displayField (VectorXd u,
                   const int N) {
  for (int i = 0; i < N; ++i)
    cout << u.segment(i * N, N).transpose() << endl; 
}

void displayField (VectorXd u,
                   const int N,
                   ofstream & stream) {
  for (int i = 0; i < N; ++i)
    stream << u.segment(i * N, N).transpose() << endl;
}
