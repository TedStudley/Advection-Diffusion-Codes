#pragma once

#include <Eigen/Dense>

using namespace Eigen;

void upwindMethod (Ref<VectorXd>,   // u
                   const int,       // N
                   const double,    // v
                   const double     // delta_t
                  );

void frommMethod (Ref<VectorXd>,    // u
                  const int,        // N
                  const double,     // v
                  const double      // delta_t
                 );

void frommVanLeer (Ref<VectorXd>,   // u
                   const int,       // N
                   const double,    // v
                   const double     // delta_t
                  );
