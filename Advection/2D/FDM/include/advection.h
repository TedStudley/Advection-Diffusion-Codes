#pragma once

#include <Eigen/Dense>

using namespace Eigen;

void upwindMethod (Ref<VectorXd>,   // u
                   const double,    // delta_t
                   const double,    // h
                   const Vector2d); // v
                   
void frommMethod (Ref<VectorXd>,    // u
                  const double,     // delta_t
                  const double,     // h
                  const Vector2d);  // v

void frommVanLeer (Ref<VectorXd>,   // u
                   const double,    // delta_t
                   const double,    // h
                   const Vector2d); // v
