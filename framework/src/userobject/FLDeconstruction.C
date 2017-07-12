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

#include "FLDeconstruction.h"

// libmesh includes
#include "libmesh/quadrature.h"

template <>
InputParameters
validParams<FLDeconstruction>()
{
  InputParameters params = validParams<SideUserObject>();
  params.addRequiredCoupledVar("variable", "The variable that will be integrated");
  params.addRequiredParam<int>("l_order", "Order of Legendre expansion");
  params.addRequiredParam<int>("f_order", "Order of Fourier expansion");
  params.addRequiredParam<std::string>("legendre_function",
    "Name of function to compute Legendre polynomial value at a point");
  params.addRequiredParam<std::string>("fourier_function",
    "Name of function to compute Fourier polynomial value at a point");
  params.addRequiredParam<int>("l_direction",
    "Direction of integration for Legendre polynomial");
  params.addRequiredParam<std::string>("aux_scalar_name",
    "Aux scalar to store the expansion coefficients");
  params.addRequiredParam<std::string>("surface_area_pp",
    "Name of post processor that calculates surface area");
  return params;
}

FLDeconstruction::FLDeconstruction(const InputParameters & parameters)
  : SideUserObject(parameters),
    MooseVariableInterface(this, false),
    _qp(0),
    _integral_value(0),
    _grad_u(coupledGradient("variable")),
    _l_order(getParam<int>("l_order")),
    _f_order(getParam<int>("f_order")),
    _legendre_function(dynamic_cast<LegendrePolynomial&>(_mci_feproblem.
      getFunction(parameters.get<std::string>("legendre_function")))),
    _fourier_function(dynamic_cast<FourierPolynomial&>(_mci_feproblem.
      getFunction(parameters.get<std::string>("fourier_function")))),
    _l_direction(getParam<int>("l_direction")),
    _aux_scalar_name(parameters.get<std::string>("aux_scalar_name")),
    _surface_area_pp(getPostprocessorValueByName(parameters.
      get<std::string>("surface_area_pp")))
{
  _num_entries = _l_order + 1;
//  _integral_value.assign(_num_entries, 0.0);
}

void
FLDeconstruction::initialize()
{
  _integral_value = 0;
}

void
FLDeconstruction::execute()
{
  _integral_value += computeIntegral();
}

Real
FLDeconstruction::getValue()
{
  gatherSum(_integral_value);
  return _integral_value;
}

void
FLDeconstruction::threadJoin(const UserObject & y)
{
  const FLDeconstruction & pps = static_cast<const FLDeconstruction &>(y);
  _integral_value += pps._integral_value;
}

Real
FLDeconstruction::computeIntegral()
{
  Real sum = 0;
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    sum += _JxW[_qp] * _coord[_qp] * computeQpIntegral();
  return sum;
}
