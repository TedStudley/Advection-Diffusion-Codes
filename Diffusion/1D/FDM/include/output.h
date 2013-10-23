#pragma once

#include <Eigen/Dense>

#include <iostream>
#include <ostream>

using namespace Eigen;
using namespace std;

void displayField (VectorXd u);

void displayField (VectorXd   u,
                   ofstream & stream);

void outputStats (const int    N,
                  const double h,
                  const double dt,
                  const int    stride,
                  const double T);

ofstream & openTeXDoc (string filename);

void makeTeXRow (VectorXd   u,
                 ofstream & stream);

void closeTeXDoc (ofstream & stream);
