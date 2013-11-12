#pragma once

#include <Eigen/Dense>

#include <iostream>
#include <ostream>

using namespace Eigen;
using namespace std;

void displayField (VectorXd);

void displayField (VectorXd,
                   ofstream &);

void outputStats (const int,
                  const double,
                  const double,
                  const int,
                  const double);

ofstream & openTeXDoc (string filename);

void makeTeXRow (VectorXd   u,
                 ofstream & stream);

void closeTeXDoc (ofstream & stream);
