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
#include "MooseVariableScalar.h"
#include "SystemBase.h"

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
  addMooseVariableDependency(mooseVariable());

  _num_entries = _l_order + 1;
  _integral_value.assign(_num_entries, 0.0);

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

void
FLDeconstruction::initialize()
{
  _integral_value.assign(_num_entries, 0.0);
}

void
FLDeconstruction::execute()
{
  for (int l = 0; l < _num_entries; ++l)
    _integral_value[l] += computeIntegral(_f_order, l);
}

Real
FLDeconstruction::getValue(int N)
{
  gatherSum(_integral_value[N]);
  return _integral_value[N];
}

void
FLDeconstruction::threadJoin(const UserObject & y)
{
  for (int l = 0; l < _num_entries; ++l)
  {
    const FLDeconstruction & pps = static_cast<const FLDeconstruction &>(y);
    _integral_value[l] += pps._integral_value[l];
  }
}

Real
FLDeconstruction::computeIntegral(int f, int l)
{
  Real sum = 0;
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    sum += _JxW[_qp] * _coord[_qp] * computeQpIntegral(f, l);
  return sum;
}

Real
FLDeconstruction::computeQpIntegral(int f, int l)
{
  Real l_func = _legendre_function.getPolynomialValue(_t,
    _q_point[_qp](_l_direction), l);
  Real f_func = _fourier_function.getPolynomialValue(_t,
    _q_point[_qp](_fdir1), _q_point[_qp](_fdir2), f);

  return _grad_u[_qp] * _normals[_qp] * l_func * f_func *
    4.0 * M_PI / _surface_area_pp;
}

void
FLDeconstruction::finalize()
{
  /* Store the result of the user object in a SCALAR variable.*/
  MooseVariableScalar & scalar =
    _fe_problem.getScalarVariable(_tid, _aux_scalar_name);
  scalar.reinit();

  std::vector<dof_id_type> & dof = scalar.dofIndices();

  for (int i = 0; i < _num_entries; ++i)
    scalar.sys().solution().set(dof[i], getValue(i));

  scalar.sys().solution().close();
}
