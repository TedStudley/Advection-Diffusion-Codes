#pragma once

#include <Eigen/Dense>

using namespace Eigen;

void squareWave (Eigen::Ref<Eigen::VectorXd>);  // u

void fourierSquare (Eigen::Ref<Eigen::VectorXd>,  // u
                    const double,                 // kappa
                    const double);                // t0

void sineWave (Eigen::Ref<Eigen::VectorXd>,
               const int);

void sineWave (Eigen::Ref<Eigen::VectorXd>,
               const int,
               const double,
               const double);
