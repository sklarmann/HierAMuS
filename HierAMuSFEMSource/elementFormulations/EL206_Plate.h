// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once


#include "elementFormulations/GenericElementFormulationInterface.h"
#include "types/MatrixTypes.h"

#include <Eigen/Dense>

namespace HierAMuS::Geometry {
struct H1Shapes;
}

namespace HierAMuS::Elementformulations {

class EL206_Plate : public GenericElementFormulationInterface<FiniteElement::Face> {
public:
  explicit EL206_Plate(PointerCollection *ptrCol);
  ~EL206_Plate() override;
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

  auto getNumberOfIntergrationPoints(PointerCollection &pointers,
                                     FiniteElement::GenericFiniteElement *elem)
      -> indexType override;

  // plot
  void toParaviewAdaper(PointerCollection &pointers,
                                FiniteElement::Face &elem,
                                vtkPlotInterface &paraviewAdapter,
                                ParaviewSwitch control) override;



private:
  Types::Matrix33<prec> getMaterialMatrix();
  Types::Matrix3X<prec>
  getLocalStrainStressInterpolation(FiniteElement::GenericFiniteElement *elem,
                                    prec xi, prec eta);
  Types::Matrix3X<prec> getBMatrix(Geometry::H1Shapes &shapesDisp, Geometry::H1Shapes &shapesRot);
  void toParaviewQuadrilateral(PointerCollection &pointers,
                               FiniteElement::GenericFiniteElement *elem,
                               vtkPlotInterface &paraviewAdapter,
                               ParaviewSwitch control);

  void setTangentResidualDispFormulation(
    PointerCollection& pointers,
    FiniteElement::Face *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  void setTangentResidualHuWashizuFormulation(
    PointerCollection& pointers,
    FiniteElement::Face *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  indexType meshIdDisp, meshIdRot, dispOrder, rotOrder;
  indexType mode;
  prec EI, GA, GI;
};

} // namespace HierAMuS
