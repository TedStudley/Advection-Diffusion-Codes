#pragma once

#include <Eigen/Dense>

using namespace Eigen;

void squareWave (Ref<VectorXd>);  // u

void fourierSquare (Ref<VectorXd>,  // u
                    const double,   // kappa
                    const double);  // t0

void sineWave (Ref<VectorXd>, // u
               const int,     // k1
               const int);    // k2

void sineWave (Ref<VectorXd>, // u
               const int,     // k1
               const int,     // k2
               const double,  // kappa
               const double); // t0
