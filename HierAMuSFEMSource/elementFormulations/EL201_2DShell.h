// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "elementFormulations/GenericElementFormulationInterface.h"


namespace HierAMuS::Geometry{
  struct H1Shapes;
}
//#include <types/MatrixTypes.h>

//#include <Eigen/Dense>


namespace HierAMuS::Elementformulations {

class EL201_2DShell : public GenericElementFormulationInterface<FiniteElement::Face> {
public:
  explicit EL201_2DShell(PointerCollection *ptrCol);
  ~EL201_2DShell() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection &pointers,
                           FiniteElement::Face &elem) override;
  void AdditionalOperations(PointerCollection& pointers, FiniteElement::Face &elem) override;
  auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> std::vector<DegreeOfFreedom *> override;
  void setTangentResidual(
    PointerCollection& pointers,
    FiniteElement::Face &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;

 
  // plot
  void toParaviewAdaper(PointerCollection &pointers,
                                FiniteElement::Face &elem,
                                vtkPlotInterface &paraviewAdapter,
                                ParaviewSwitch control) override;

  auto getHistoryDataStructure() -> const HistoryDataStructure & override;
  auto getNumberOfIntergrationPoints(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> indexType override;
  
private:
  Types::Matrix3X<prec>
  getLocalStrainStressInterpolation(PointerCollection& pointers,
                                    FiniteElement::Face &elem, prec xi, prec eta);
  static Types::Matrix3X<prec> getBMatrix(const Types::Matrix2X<prec> &shapeDeriv,
                                   indexType numDofs);
  static Types::VectorXT<prec> getHMatrix(const Types::VectorX<prec>& shapes,
      indexType numDofs);
  static Types::VectorXT<prec> getBvolMatrix(const Types::Matrix2X<prec>& shapeDeriv,
      indexType numDofs);
  // void toParaviewQuadrilateral(PointerCollection &pointers,
  //                              FiniteElement::GenericFiniteElement *elem,
  //                              vtkPlotInterface &paraviewAdapter,
  //                              ParaviewSwitch control);

  void setTangentResidualDispFormulation(
    PointerCollection& pointers,
    FiniteElement::Face &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  void setTangentResidualDispNonLinear(
    PointerCollection& pointers,
    FiniteElement::Face &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  void setTangentResidualHuWashizuFormulation(
    PointerCollection& pointers,
    FiniteElement::Face &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  void setTangentResidualDispAndStreFormulation(
    PointerCollection& pointers,
    FiniteElement::Face &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic>& stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1>& residual, std::vector<DegreeOfFreedom*>& Dofs);

  auto getDeformationGradient(Geometry::H1Shapes &shapes, Types::VectorX<prec> &solution) -> Types::Matrix22<prec>;

  auto getNonLinearBMatrix(Geometry::H1Shapes &shapes, Types::Matrix22<prec> &F) -> Types::Matrix3X<prec>;

  auto getGeometricMatrix(Geometry::H1Shapes &shapes, Types::VectorX<prec> &stress) -> Types::MatrixXX<prec>;

  auto getIntegrationPoints(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> IntegrationPoints;

  indexType meshIdDisp, intOrderDisp;
  indexType meshIdStre, intOrderStre;
  indexType mode;


  const static HistoryDataStructure m_HistoryDataStructure;
};

} // namespace HierAMuS
