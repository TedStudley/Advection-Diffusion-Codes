#pragma once

#include <Eigen/Dense>

#include <iostream>
#include <ostream>

using namespace Eigen;
using namespace std;

void displayField (VectorXd,
                   const int);

void displayField (VectorXd,
                   const int,
                   ostream &);
