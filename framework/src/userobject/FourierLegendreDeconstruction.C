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

#include "FourierLegendreDeconstruction.h"
#include "MooseVariableScalar.h"
#include "SystemBase.h"

template<>
InputParameters validParams<FourierLegendreDeconstruction>()
{
  InputParameters params = validParams<SideIntegralUserObject>();
  params.addRequiredCoupledVar("variable",
    "The name of the variable that will be integrated");
  params.addRequiredParam<std::string>("legendre_function",
    "Name of function to compute the Legendre polynomial at a point");
  params.addRequiredParam<std::string>("fourier_function",
    "Name of function to compute the Fourier polynomial at a point");
  params.addRequiredParam<int>("l_direction",
    "Direction of integration for Legendre polynomial");
  params.addRequiredParam<int>("l_order", "The order of the Legendre expansion");
  params.addRequiredParam<int>("f_order", "The order of the Fourier expansion");
  params.addRequiredParam<std::string>("aux_scalar_name",
    "Aux scalar to store the Legendre expansion coefficient");
  params.addRequiredParam<std::string>("surface_area_pp",
    "The name of the post processor that calculates surface area");
  return params;
}

FourierLegendreDeconstruction::FourierLegendreDeconstruction(const InputParameters & parameters) :
    SideIntegralUserObject(parameters),
    MooseVariableInterface(this, false),
    _u(coupledValue("variable")),
    _legendre_function(dynamic_cast<LegendrePolynomial&>(_mci_feproblem.
      getFunction(parameters.get<std::string>("legendre_function")))),
    _fourier_function(dynamic_cast<FourierPolynomial&>(_mci_feproblem.
      getFunction(parameters.get<std::string>("fourier_function")))),
    _l_direction(parameters.get<int>("l_direction")),
    _l_order(parameters.get<int>("l_order")),
    _f_order(parameters.get<int>("f_order")),
    _aux_scalar_name(parameters.get<std::string>("aux_scalar_name")),
    _surface_area_pp(getPostprocessorValueByName(parameters.
      get<std::string>("surface_area_pp")))

{
  addMooseVariableDependency(mooseVariable());

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
}

Real
FourierLegendreDeconstruction::computeQpIntegral()
{
  Real l_func = _legendre_function.
    getPolynomialValue(_t, _q_point[_qp](_l_direction), _l_order);
  Real f_func = _fourier_function.
    getPolynomialValue(_t, _q_point[_qp](_fdir1), _q_point[_qp](_fdir2), _f_order);
  return _u[_qp] * l_func * f_func * 4.0 * M_PI / _surface_area_pp;
}
void
FourierLegendreDeconstruction::finalize()
{
  MooseVariableScalar & scalar =
    _fe_problem.getScalarVariable(_tid, _aux_scalar_name);
  scalar.reinit();
  std::vector<dof_id_type> & dof = scalar.dofIndices();
  scalar.sys().solution().set(dof[_l_order], getValue());
  scalar.sys().solution().close();
}
