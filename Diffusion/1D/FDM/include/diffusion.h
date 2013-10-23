#pragma once

#include <Eigen/Dense>

using namespace Eigen;

void forwardEuler (Ref<VectorXd> u,   
                   const double  dt,    
                   const double  h,    
                   const double  kappa);

void backwardEuler (Ref<VectorXd> u,  
                    const double  dt,   
                    const double  h,   
                    const double  kappa);

void crankNicolson (Ref<VectorXd> u,
                    const double  dt,
                    const double  h,
                    const double  kappa);

void BDF2 (Ref<VectorXd> u,
           Ref<VectorXd> u1,
           const double  dt,
           const double  h,
           const double  kappa);

void TR_BDF2 (Ref<VectorXd> u,  
              const double  dt, 
              const double  h,  
              const double  kappa); 
