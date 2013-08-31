#pragma once

#include <Eigen/Dense>

using namespace Eigen;

void forwardEuler (Ref<VectorXd>,   // u
                   const int,       // N
                   const double,    // kappa
                   const double     // delta_t
                  );

void backwardEuler (Ref<VectorXd>,  // u
                    const int,      // N
                    const double,   // kappa
                    const double    // delta_t
                   );

void crankNicolson (Ref<VectorXd>,  // u
                     const int,     // N
                     const double,  // kappa
                     const double,  // delta_t
                     const double c = 0.5
                    );

void BDF2 (Ref<VectorXd>, // u
           Ref<VectorXd>, // u1
           const int,     // N
           const double,  // kappa
           const double); // delta_t
