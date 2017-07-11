/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "FourierPolynomial.h"

template<>
InputParameters validParams<FourierPolynomial>()
{
  InputParameters params = validParams<Function>();
  params.addRequiredParam<std::vector<Real>>("center",
    "center coordinates of circle for normalization.");
  params.addParam<bool>("dbg", "Print debug output");
  return params;
}

FourierPolynomial::FourierPolynomial(const InputParameters & parameters) :
    Function(parameters),
    _center(parameters.get<std::vector<Real>>("center")),
    _dbg(parameters.get<bool>("dbg"))
{
}

FourierPolynomial::~FourierPolynomial()
{
}

Real
FourierPolynomial::value(Real t, const Point & p)
{
  mooseWarning("value() in FourierPolynomial should not be used");
  return 0.0;
}

Real
FourierPolynomial::getPolynomialValue(Real t, Real p1, Real p2, int n)
{
  Real y = p2 - _center[1];
  Real x = p1 - _center[0];

  Real phi = atan2(p2, p1);
  if (n == 0)
    return 1.0 / sqrt(2.0 * M_PI);
  else
    return cos(n * phi) / sqrt(M_PI);
}
