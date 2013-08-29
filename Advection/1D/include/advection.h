#pragma once

#include <Eigen/Dense>

using namespace Eigen;

void upwindMethod (Ref<VectorXd>,   // u
                   unsigned int,    // N
                   double,          // v
                   double           // delta_t
                  );

void frommMethod (Ref<VectorXd>,    // u
                  unsigned int,     // N
                  double,           // v
                  double            // delta_t
                 );

void frommVanLeer (Ref<VectorXd>,   // u
                   unsigned int,    // N
                   double,          // v
                   double           // delta_t
                  );
