// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "geometry/Base.h"
#include "geometry/GeometryTypes.h"
#include <forwarddeclaration.h>


#include <elementFormulations/GenericElementFormulation.h>
#include <types/MatrixTypes.h>

#include <Eigen/Dense>

namespace HierAMuS::Elementformulations {

class EL206_Plate : public GenericElementFormulation {
public:
  explicit EL206_Plate(PointerCollection *ptrCol);
  ~EL206_Plate() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) override;
  void AdditionalOperations(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) override;
  auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> std::vector<DegreeOfFreedom *> override;
  void setTangentResidual(
    PointerCollection& pointers,
    FiniteElement::GenericFiniteElement *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;

  auto getNumberOfIntergrationPoints(PointerCollection &pointers,
                                     FiniteElement::GenericFiniteElement *elem)
      -> indexType override;

  // plot
  void toParaviewAdaper(PointerCollection &pointers,
                                FiniteElement::GenericFiniteElement *elem,
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
    FiniteElement::GenericFiniteElement *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  void setTangentResidualHuWashizuFormulation(
    PointerCollection& pointers,
    FiniteElement::GenericFiniteElement *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  indexType meshIdDisp, meshIdRot, dispOrder, rotOrder;
  indexType mode;
  prec EI, GA, GI;
};

} // namespace HierAMuS
