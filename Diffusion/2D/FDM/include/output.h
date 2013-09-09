#pragma once

#include <Eigen/Dense>

#include <fstream>

using namespace Eigen;
using namespace std;


void displayField (VectorXd);

void displayField (VectorXd,
                   ofstream &);
