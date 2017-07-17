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

#ifndef FOURIERLEGENDRERECONSTRUCTION_H
#define FOURIERLEGENDRERECONSTRUCTION_H

#include "Function.h"
#include "math.h"
#include "FourierPolynomial.h"
#include "LegendrePolynomial.h"

class FourierLegendreReconstruction : public Function
{
public:
  FourierLegendreReconstruction(const InputParameters & parameters);
  virtual ~FourierLegendreReconstruction();

  virtual Real value(Real t, const Point & p);

protected:

  bool _dbg;
  int _l_direction;
  int _fdir1;
  int _fdir2;
  int _l_order;
  int _f_order;
  std::vector<VariableValue *> _poly_coeffs;
  LegendrePolynomial & _legendre_function;
  FourierPolynomial & _fourier_function;

private:

};

template<>
InputParameters validParams<FourierLegendreReconstruction>();

#endif
