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

#ifndef FLDECONSTRUCTION_H
#define FLDECONSTRUCTION_H

// MOOSE includes
#include "SideUserObject.h"
#include "MooseVariableInterface.h"
#include "LegendrePolynomial.h"
#include "FourierPolynomial.h"

class FLDeconstruction;

template <>
InputParameters validParams<FLDeconstruction>();

class FLDeconstruction : public SideUserObject, public MooseVariableInterface
{
public:
  FLDeconstruction(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

  virtual Real getValue(int);

protected:
  virtual Real computeQpIntegral(int, int);
  virtual Real computeIntegral(int, int);

  unsigned int _qp;
  const VariableGradient & _grad_u;
  int _l_order;
  int _f_order;
  int _num_entries;
  LegendrePolynomial & _legendre_function;
  FourierPolynomial & _fourier_function;
  int _l_direction;
  int _fdir1;
  int _fdir2;

  std::string _aux_scalar_name;
  const PostprocessorValue & _surface_area_pp;

  std::vector<Real> _integral_value;
};

#endif
