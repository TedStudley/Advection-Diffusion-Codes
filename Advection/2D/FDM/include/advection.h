#pragma once

#include <Eigen/SparseCore>

using namespace Eigen;

void upwindMethod (Ref<VectorXd>  u, 
                   const double   dt,
                   const double   h,
                   const Vector2d v);
                   
void frommMethod (Ref<VectorXd>  u,
                  const double   dt,
                  const double   h, 
                  const Vector2d v);

void frommVanLeer (Ref<VectorXd>  u,
                   const double   dt,
                   const double   h, 
                   const Vector2d v);

