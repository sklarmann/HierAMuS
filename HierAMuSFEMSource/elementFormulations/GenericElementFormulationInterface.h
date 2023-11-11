// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "GenericElementFormulation.h"

class managementClass;
namespace HierAMuS {
class PointerCollection;
class ParameterList;
class DegreeOfFreedom;
class vtkPlotInterface;
enum class ParaviewSwitch;
class ParameterList;
namespace FiniteElement {
class GenericFiniteElement;
class Edge;
class Face;
class Volume;
class beamInterfaceElement3D;
} // namespace FiniteElement

namespace Elementformulations {

template <class T>
class GenericElementFormulationInterface : public GenericElementFormulation {
public:
  explicit GenericElementFormulationInterface(PointerCollection *pointer)
      : GenericElementFormulation(pointer){};

  
  virtual void setDegreesOfFreedom(PointerCollection &pointers,
                                   T &elem) = 0;

  virtual void AdditionalOperations(PointerCollection &pointers,
                                    T &elem) = 0;

  virtual void setTangentResidual(PointerCollection &pointers, T &elem,
                                  Types::MatrixXX<prec> &stiffness,
                                  Types::VectorX<prec> &residual,
                                  std::vector<DegreeOfFreedom *> &Dofs) = 0;

  virtual void toParaviewAdaper(PointerCollection &pointers, T &elem,
                                vtkPlotInterface &paraviewAdapter,
                                ParaviewSwitch control) = 0;
};
} // namespace Elementformulations
} // namespace HierAMuS
