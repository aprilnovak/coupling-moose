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

#include "FourierLegendreReconstruction.h"
#include <iostream>

template<>
InputParameters validParams<FourierLegendreReconstruction>()
{
  InputParameters params = validParams<Function>();
  params.addParam<bool>("dbg", false, "Print debug output");
  params.addRequiredParam<int>("l_order", "The order of the Legendre expansion");
  params.addRequiredParam<int>("f_order", "The order of the Fourier expansion");
  params.addRequiredParam<int>("l_direction", "The coordinate direction in which"
    " the Legendre expansion is applied (0 = x, 1 = y, 2 = z).");
  params.addRequiredParam<std::string>("legendre_function", "Name of the function"
    " to compute the Legendre polynomial at a point");
  params.addRequiredParam<std::string>("fourier_function", "Name of the function"
    " to compute the Fourier polynomial at a point");
  params.addRequiredCoupledVar("poly_coeffs", "Name of the aux variables"
    " containing the Zernike coefficients");
  return params;
}

FourierLegendreReconstruction::FourierLegendreReconstruction(const InputParameters & parameters) :
    Function(parameters),
    _dbg(parameters.get<bool>("dbg")),
    _l_direction(parameters.get<int>("l_direction")),
    _l_order(parameters.get<int>("l_order")),
    _f_order(parameters.get<int>("f_order")),
    _legendre_function(dynamic_cast<LegendrePolynomial&>(_mci_feproblem.
      getFunction(parameters.get<std::string>("legendre_function")))),
    _fourier_function(dynamic_cast<FourierPolynomial&>(_mci_feproblem.
      getFunction(parameters.get<std::string>("fourier_function"))))
{
  if (_l_direction == 0) // Legendre in x-direction, Zernike in y-z
  {
    _fdir1 = 1;
    _fdir2 = 2;
  }
  else if (_l_direction == 1) // Legendre in y-direction, Zernike in x-z
  {
    _fdir1 = 0;
    _fdir2 = 2;
  }
  else // Legendre in z-direction, Zernike in x-y
  {
    _fdir1 = 0;
    _fdir2 = 1;
  }

  // Get the coupled scalar variables storing the expansion coefficients. This
  // will give an error if the f_order is too large (not enough coupled
  // variables).
  // TODO: implement check that each coupled variable is the same size and
  //       that the number of coefficients per variable match expected
  //       for Legendre order
  for (int i = 0; i <= _f_order; i++)
  {
    _poly_coeffs.push_back(&coupledScalarValue("poly_coeffs", i));
  }

  // TODO: check to see if the number of provided variables exceeds the number
  // expected based on the provided Fourier order (i.e. if f_order was made
  // too small).
}

FourierLegendreReconstruction::~FourierLegendreReconstruction()
{
}

Real
FourierLegendreReconstruction::value(Real t, const Point & p)
{
  Real val = 0.0;

  for (int f = 0; f <= _f_order; ++f)
  {
    for (int l = 0; l <= _l_order; ++l)
    {
      val += (*_poly_coeffs[f])[l]
        * _legendre_function.getPolynomialValue(t, p(_l_direction), l)
        * _fourier_function.getPolynomialValue(t, p(_fdir1), p(_fdir2), f);
    }
  }

  return val;
}
