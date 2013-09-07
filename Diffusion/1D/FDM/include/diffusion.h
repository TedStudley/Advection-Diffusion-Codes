#pragma once

#include <Eigen/Dense>

using namespace Eigen;

void forwardEuler (Ref<VectorXd>,   // u
                   const int,       // N
                   const double);   // mu

void backwardEuler (Ref<VectorXd>,  // u
                    const int,      // N
                    const double);  // mu

void crankNicolson (Ref<VectorXd>,  // u
                     const int,     // N
                     const double); // mu

void BDF2 (Ref<VectorXd>, // u
           Ref<VectorXd>, // u1
           const int,     // N
           const double); // mu

void TR_BDF2 (Ref<VectorXd>,  // u
              const int,      // N
              const double,   // delta_t
              const double);  // mu
