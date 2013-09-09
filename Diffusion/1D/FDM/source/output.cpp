#include <output.h>

#include <Eigen/Dense>

#include <iostream>
#include <fstream>

using namespace Eigen;
using namespace std;


void displayField (VectorXd u) {
  cout << u << endl;
}

void displayField (VectorXd u,
                   ofstream & stream) {
  stream << u << endl;
}  
