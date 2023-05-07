// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <forwarddeclaration.h>

#include <elementFormulations/GenericElementFormulation.h>
#include <Eigen/Dense>
#include <types/MatrixTypes.h>
#include <vector>

namespace HierAMuS::Elementformulations {

class EL204_BeamInterface2D : public GenericElementFormulation {
public:
  explicit EL204_BeamInterface2D(PointerCollection *ptrCol);
  ~EL204_BeamInterface2D() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) override;
  void AdditionalOperations(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) override;
  auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> std::vector<DegreeOfFreedom *> override;
  void setTangentResidual(PointerCollection& pointers,
                          FiniteElement::GenericFiniteElement *elem,
                          Types::MatrixXX<prec> &stiffness,
                          Types::VectorX<prec> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;

  void setTangentResidualMixed(FiniteElement::GenericFiniteElement *elem,
                               Types::MatrixXX<prec> &stiffness,
                               Types::VectorX<prec> &residual,
                               std::vector<DegreeOfFreedom *> &Dofs);

  auto getNumberOfIntergrationPoints(PointerCollection &pointers,
                                     FiniteElement::GenericFiniteElement *elem)
      -> indexType override;

  //Plot
  void toParaviewAdaper(PointerCollection &pointers,
                        FiniteElement::GenericFiniteElement *elem,
                        vtkPlotInterface &paraviewAdapter,
                        ParaviewSwitch control) override;

private:
  indexType meshIdDisp, intOrder, planeStrain;
  indexType meshIdWarp, meshIdRot;
  prec emod, nu, thick;
  indexType mode;
};

} // namespace HierAMuS
