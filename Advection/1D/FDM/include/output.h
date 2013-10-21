#pragma once

#include <Eigen/Dense>

#include <iostream>
#include <ostream>

using namespace Eigen;
using namespace std;

void displayField (VectorXd);

void displayField (VectorXd,
                   ostream &);

void outputStats (const int,
                  const double,
                  const double,
                  const int,
                  const double);

ofstream & openTeXDoc (string);

void makeTeXRow (ofstream &,
                 VectorXd);

void closeTeXDoc (ofstream &);
