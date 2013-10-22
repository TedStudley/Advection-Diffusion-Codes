#pragma once

#include <Eigen/Dense>

using namespace Eigen;

void squareWave    (Ref<VectorXd> u);  // u

void fourierSquare (Ref<VectorXd> u,  // u
                    const double  kappa,   // kappa
                    const double  t0);  // t0

void sineWave      (Ref<VectorXd> u, // u
               const int);    // k

void sineWave (Ref<VectorXd>, // u
               const int,     // k
               const double,  // kappa
               const double); // t0
