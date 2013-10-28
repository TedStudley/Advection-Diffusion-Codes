#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;

template<typename Scalar> struct buildDelta {
  EIGEN_EMPTY_STRUCT_CTOR(buildDelta)
  typedef double result_type;
  double operator() (const Scalar& a, const Scalar& b) const { return ((b > 0) ? a : 0); }
};

template <typename Scalar> struct unarySign {
  EIGEN_EMPTY_STRUCT_CTOR(unarySign)
  typedef double result_type;
  double operator() (const Scalar& a) const { return ((a < 0) ? -1 : (a == 0) ? 0 : 1); }
};

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
