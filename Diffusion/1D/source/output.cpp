#include <output.h>

#include <Eigen/Dense>

#include <iostream>
#include <fstream>

using namespace Eigen;
using namespace std;


void displayField (VectorXd u,
                   const int N) {
  cout << u.transpose () << endl;
}

void displayField (VectorXd u,
                   const int N,
                   ofstream & stream) {
  stream << u.transpose() << endl;
}  
