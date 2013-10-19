#pragma once

#include <Eigen/Dense>

#include <fstream>

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
