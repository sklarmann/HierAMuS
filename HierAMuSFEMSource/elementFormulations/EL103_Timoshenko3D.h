// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "types/MatrixTypes.h"

#include "elementFormulations/GenericElementFormulationInterface.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <Eigen/Dense>

namespace HierAMuS {
namespace FiniteElement {
class Edge;
}
namespace Elementformulations {

class EL103_Timoshenko3D : public GenericElementFormulationInterface<FiniteElement::Edge> {
public:
  explicit EL103_Timoshenko3D(PointerCollection *ptrCol);
  ~EL103_Timoshenko3D() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection &pointers,
                           FiniteElement::Edge &elem) override;
  void AdditionalOperations(PointerCollection &pointers,
                            FiniteElement::Edge &elem) override;
  void updateRVEHistory(PointerCollection &pointers,
                        FiniteElement::GenericFiniteElement *elem) override;
  auto getDofs(PointerCollection &pointers,
               FiniteElement::GenericFiniteElement *elem)
      -> std::vector<DegreeOfFreedom *> override;
  void setTangentResidual(
      PointerCollection &pointers, FiniteElement::Edge &elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs) override;
  void setMass(PointerCollection &pointers,
               FiniteElement::GenericFiniteElement *elem,
               Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
               Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
               std::vector<DegreeOfFreedom *> &Dofs) override;
  void getElementsLocalNodalReactions(
      PointerCollection &ptrCol, FiniteElement::GenericFiniteElement *elem,
      std::map<indexType, std::vector<prec>> &vReacs) override;

  // plot
  void toParaviewAdaper(PointerCollection &pointers, FiniteElement::Edge &elem,
                        vtkPlotInterface &paraviewAdapter,
                        ParaviewSwitch control) override;

  auto getHistoryDataStructure() -> const HistoryDataStructure & override;
  auto getNumberOfIntergrationPoints(PointerCollection &pointers,
                                     FiniteElement::GenericFiniteElement *elem)
      -> indexType override;

private:
  void setTangentResidualLinear(
      PointerCollection &pointers, FiniteElement::Edge &elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);
  void setTangentResidualLinearFE2(
      PointerCollection &pointers, FiniteElement::Edge &elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);
  void setTangentResidualLinearThermal(
      PointerCollection &pointers, FiniteElement::Edge &elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);
  void setTangentResidualLinearThermalFE2(
      PointerCollection &pointers, FiniteElement::Edge &elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);
  void setTangentResidualNonlinear(
      PointerCollection &pointers, FiniteElement::Edge &elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);

  auto getLinearCMat() -> Types::Matrix66<prec>;
  auto getLinearCMatThermal() -> Types::MatrixXX<prec>;

  auto getIntegrationPoints(PointerCollection &pointers,
                            FiniteElement::GenericFiniteElement *elem)
      -> IntegrationPoints;

  Types::Matrix33<prec> getMaterialMatrix() const;
  indexType meshIdDisp, meshIdRot, meshIdTemp, disporder, rotorder, temporder;
  indexType geo, thermal, mode, fe2, RVECmat;
  prec E, G, A, Ix, Iy, Iz, ky, kz, alpha, kappa;
  bool m_RVECmatInit;
  Types::MatrixXX<prec> m_CMat;

  const static HistoryDataStructure m_HistoryDataStructure;
};

} // namespace Elementformulations
} // namespace HierAMuS
