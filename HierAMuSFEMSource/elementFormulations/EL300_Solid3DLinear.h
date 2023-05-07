// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once

#include "MatrixTypes.h"
#include "geometry/Base.h"
#include "geometry/GeometryTypes.h"
#include <forwarddeclaration.h>

#include <elementFormulations/GenericElementFormulation.h>
#include <Eigen/Dense>

namespace HierAMuS {
namespace Elementformulations {

class EL300_Solid3DLinear : public GenericElementFormulation {
public:
  EL300_Solid3DLinear(PointerCollection *ptrCol);
  ~EL300_Solid3DLinear() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) override;
  void AdditionalOperations(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) override;
  auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> std::vector<DegreeOfFreedom*> override;
  void setTangentResidual(
    PointerCollection& pointers,
    FiniteElement::GenericFiniteElement *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;

  // Paraview
  void toParaviewAdaper(PointerCollection &pointers,
                                FiniteElement::GenericFiniteElement *elem,
                                vtkPlotInterface &paraviewAdapter,
                                ParaviewSwitch control) override;
  auto getHistoryDataStructure() -> const HistoryDataStructure & override;
  auto getNumberOfIntergrationPoints(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> indexType override;

private:
  void setTangentResidualLinear(
    PointerCollection& pointers,
    FiniteElement::GenericFiniteElement *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);
  void setTangentResidualNonLinear(
    PointerCollection& pointers,
    FiniteElement::GenericFiniteElement *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  auto getDeformationGradient(Geometry::H1Shapes &shapes, Types::VectorX<prec> &solution) -> Types::Matrix33<prec>;
  auto getNonLinearBMatrix(Geometry::H1Shapes &shapes, Types::Matrix33<prec> &F) -> Types::Matrix6X<prec>;
  auto getLinearBMatrix(Geometry::H1Shapes &shapes) -> Types::Matrix6X<prec>;
  auto getGeometricMatrix(Geometry::H1Shapes &shapes, Types::VectorX<prec> &stress) -> Types::MatrixXX<prec>;

  auto getBTCBLinear(Geometry::H1Shapes &shapes,Types::MatrixXX<prec> &C, indexType nodeI, indexType nodeJ) -> Types::Matrix33<prec>;

  indexType meshIdDisp, intOrderDisp, mode;

  const static HistoryDataStructure m_HistoryDataStructure;
  //prec mu1, mu2;
  //prec e31, e33, e15;
};

} // namespace Elementformulations
} // namespace HierAMuS
