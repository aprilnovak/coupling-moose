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

#ifndef FOURIERLEGENDREDECONSTRUCTION_H
#define FOURIERLEGENDREDECONSTRUCTION_H

#include "SideIntegralUserObject.h"
#include "MooseVariableInterface.h"
#include "LegendrePolynomial.h"
#include "FourierPolynomial.h"

class FourierLegendreDeconstruction;

template<>
InputParameters validParams<FourierLegendreDeconstruction>();

class FourierLegendreDeconstruction : public SideIntegralUserObject,
  public MooseVariableInterface
{
public:
  FourierLegendreDeconstruction(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;
  virtual void finalize() override;

  const VariableValue & _u;
  LegendrePolynomial & _legendre_function;
  FourierPolynomial & _fourier_function;
  int _l_direction;
  int _fdir1;
  int _fdir2;
  int _l_order;
  int _f_order;
  std::string _aux_scalar_name;
  const PostprocessorValue & _surface_area_pp;
};

#endif
