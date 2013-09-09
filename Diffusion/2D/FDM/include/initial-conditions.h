#pragma once

#include <Eigen/Dense>

using namespace Eigen;

void squareWave (Ref<VectorXd>);  // u

void fourierSquare (Ref<VectorXd>,  // u
                    const double,   // kappa
                    const double);  // t0
