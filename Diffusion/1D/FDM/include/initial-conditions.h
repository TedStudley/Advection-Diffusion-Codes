#pragma once

#include <Eigen/Dense>

using namespace Eigen;

void squareWave (Ref<VectorXd>,   // u
                 const int);      // N

void fourierSquare (Ref<VectorXd>,  // u
                    const int,      // N
                    const double,   // kappa
                    const double);  // t0
