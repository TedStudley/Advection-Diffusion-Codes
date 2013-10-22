#pragma once

#include <Eigen/Dense>

using namespace Eigen;

void upwindMethod (Ref<VectorXd>  u,   
                   const double   dt,   
                   const double   dx,   
                   const VectorXd v); 

void frommMethod  (Ref<VectorXd>  u,   
                   const double   dt,   
                   const double   dx,   
                   const VectorXd v); 

void beamWarming  (Ref<VectorXd>  u,
                   const double   dt,
                   const double   dx,
                   const VectorXd v);

void laxWendroff  (Ref<VectorXd>  u,
                   const double   dt,
                   const double   dx,
                   const VectorXd v);

void frommVanLeer (Ref<VectorXd>  u,  
                   const double   dt,
                   const double   dx,
                   const VectorXd v);
