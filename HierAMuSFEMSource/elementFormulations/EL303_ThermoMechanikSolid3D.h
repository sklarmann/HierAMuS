// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once

#include "elementFormulations/GenericElementFormulationInterface.h"

namespace HierAMuS {
namespace Geometry {
struct H1Shapes;
}
namespace Elementformulations {

class EL303_ThermoMechanikSolid3D : public GenericElementFormulationInterface<FiniteElement::Volume> {
public:
  EL303_ThermoMechanikSolid3D(PointerCollection *ptrCol);
  ~EL303_ThermoMechanikSolid3D() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::Volume &elem) override;
  void AdditionalOperations(PointerCollection& pointers, FiniteElement::Volume &elem) override;
  auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> std::vector<DegreeOfFreedom*> override;
  void setTangentResidual(
    PointerCollection& pointers,
    FiniteElement::Volume &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;

  // Paraview
  void toParaviewAdaper(PointerCollection &pointers,
                                FiniteElement::Volume &elem,
                                vtkPlotInterface &paraviewAdapter,
                                ParaviewSwitch control) override;
  auto getHistoryDataStructure() -> const HistoryDataStructure & override;
  auto getNumberOfIntergrationPoints(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> indexType override;

private:
  void setTangentResidualLinear(
    PointerCollection& pointers,
    FiniteElement::Volume &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);
  void setTangentResidualNonLinear(
    PointerCollection& pointers,
    FiniteElement::Volume &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  auto getDeformationGradient(Geometry::H1Shapes &shapes, Types::VectorX<prec> &solution) -> Types::Matrix33<prec>;
  auto getNonLinearBMatrix(Geometry::H1Shapes &shapes, Types::Matrix33<prec> &F) -> Types::Matrix6X<prec>;
  auto getLinearBMatrix(Geometry::H1Shapes &shapes) -> Types::Matrix6X<prec>;
  auto getLinearBMatrixTrace(Geometry::H1Shapes &shapes) -> Types::VectorX<prec>;
  auto getGeometricMatrix(Geometry::H1Shapes &shapes,
                          Types::VectorX<prec> &stress)
      -> Types::MatrixXX<prec>;

  auto getBTCBLinear(Geometry::H1Shapes &shapes,Types::MatrixXX<prec> &C, indexType nodeI, indexType nodeJ) -> Types::Matrix33<prec>;

  indexType meshIdDisp, meshidTemp, intOrderDisp, mode;

  const static HistoryDataStructure m_HistoryDataStructure;
  //prec mu1, mu2;
  //prec e31, e33, e15;
  prec mu, lambda;
  prec alpha, c, rho0, T0, kappa;
};

} // namespace Elementformulations
} // namespace HierAMuS
