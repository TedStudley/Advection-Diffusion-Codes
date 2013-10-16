#pragma once

// @sect3{Include files}
#include <deal.II/base/function.h>

#include <cmath>


using namespace dealii;

// @sect3{InitialConditions namespace}
namespace InitialConditions {
  // @sect4{Square Wave Initial Condition}
  template <int dim>
  class SquareWave : public Function<dim>
  {
    public:
      SquareWave () : Function<dim>() {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int  component = 0) const;
  };

  template <int dim>
  double SquareWave<dim>::value (const Point<dim> &p,
                                 const unsigned int /* component */) const
  {
    double return_value = 1;
    for (unsigned int d = 0; d < dim; ++d)
      return_value *= ((0.25 <= p[d] && p[d] <= 0.75) ? 1 : 0);

    return return_value;
  }

  
  // @sect4{Sine Wave Initial Condition}
  template <int dim>
  class SineWave : public Function<dim>
  {
    public:
      SineWave () : Function<dim>() {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int  component = 0) const;
  };

  template <int dim>
  double SineWave<dim>::value (const Point<dim> &p,
                               const unsigned int /* component */) const
  {
    double return_value = 1;
    for (unsigned int d = 0; d < dim; ++d)
      return_value *= std::sin (2 * M_PI * p[d]);

    return return_value;
  }
}
