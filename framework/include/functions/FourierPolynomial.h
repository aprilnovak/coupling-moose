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

#ifndef FOURIERPOLYNOMIAL_H
#define FOURIERPOLYNOMIAL_H

#include "Function.h"
#include "math.h"

class FourierPolynomial : public Function
{
public:
  FourierPolynomial(const InputParameters & parameters);
  virtual ~FourierPolynomial();

  virtual Real value(Real t, const Point & p);
  Real getPolynomialValue(Real t, Real p1, Real p2,  int n);

protected:

private:
  std::vector<Real> _center;
  bool _dbg;                     // Debug flag to print debug output
};

template<>
InputParameters validParams<FourierPolynomial>();

#endif
