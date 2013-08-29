#pragma once

#include <Eigen/Dense>

using namespace Eigen;

void squareWave (Ref<VectorXd>,   // u
                 unsigned int);   // N

void fourierSquare (Ref<VectorXd>,  // u
                    unsigned int,   // N
                    double);        // t0
