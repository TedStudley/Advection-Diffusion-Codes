#include <output.h>

#include <Eigen/Dense>

#include <iostream>
#include <ostream>

using namespace Eigen;
using namespace std;


void displayField (VectorXd u,
                   const int N) {
  cout << u.transpose() << endl;
}

void displayField (VectorXd u,
                   const int N,
                   ostream & stream) {
  stream << u.transpose() << endl;
}
