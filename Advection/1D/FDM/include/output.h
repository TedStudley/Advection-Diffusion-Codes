#pragma once

#include <Eigen/Dense>

#include <iostream>
#include <ostream>

using namespace Eigen;
using namespace std;

void displayField (VectorXd u);

void displayField (VectorXd  u,
                   ostream & stream);

void outputStats (const int N,
                  const double,
                  const double,
                  const int,
                  const double);

ofstream & openTeXDoc (string filename);

void makeTeXRow (VectorXd u,
                 ofstream & stream);

void closeTeXDoc (ofstream & stream);
