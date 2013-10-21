#pragma once

#include <Eigen/Dense>

using namespace Eigen;

void upwindMethod (Ref<VectorXd>,   // u
                   const double,    // delta_t
                   const double,    // h
                   const VectorXd); // v

void frommMethod (Ref<VectorXd>,    // u
                  const double,     // delta_t
                  const double,     // h
                  const VectorXd);  // v

void frommVanLeer (Ref<VectorXd>,   // u
                   const double,    // delta_t
                   const double,    // h
                   const VectorXd); // v
