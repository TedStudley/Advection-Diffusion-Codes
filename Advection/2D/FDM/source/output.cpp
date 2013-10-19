#include <output.h>

#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;
using namespace std;

void displayField (VectorXd u) {
  const int N = sqrt (u.rows ());

  for (int i = 0; i < N; ++i)
    cout << u.segment(i * N, N).transpose() << endl; 
}

void displayField (VectorXd u,
                   ofstream & stream) {
  const int N = sqrt (u.rows ());

  for (int i = 0; i < N; ++i)
    stream << u.segment(i * N, N).transpose() << endl;
}

void outputStats (const int N,
                  const double h,
                  const double dt,
                  const int stride,
                  const double T) {
  std::cerr << "==> Run statistics:" << endl
            << "====> Subdivisions (N)       = " << N << endl
            << "====> Timestep (dt)          = " << dt << endl
            << "====> Output stride (stride) = " << stride << endl
            << "====> End time (T)           = " << T << endl << endl;
}

