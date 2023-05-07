// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "MatrixTypes.h"
#include <forwarddeclaration.h>

#include <elementFormulations/GenericElementFormulation.h>
#include <Eigen/Dense>

#include "control/ParameterList.h"

#include "datatypes.h"
#include "finiteElements/NormTypes.h"


namespace HierAMuS::Elementformulations {


  class EL205_HDivTest : public GenericElementFormulation {
    public:
      explicit EL205_HDivTest(PointerCollection *ptrCol);
      ~EL205_HDivTest() override;
      void readData(PointerCollection &pointers, ParameterList &list) override;
      void setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) override;
      void AdditionalOperations(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) override;
      auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
      -> std::vector<DegreeOfFreedom *> override;

      void setTangentResidual(PointerCollection& pointers,
                              FiniteElement::GenericFiniteElement* elem,
                              Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic>& stiffness,
                              Eigen::Matrix<prec, Eigen::Dynamic, 1>& residual, std::vector<DegreeOfFreedom*>& Dofs) override;

      auto
      getNumberOfIntergrationPoints(PointerCollection &pointers,
                                    FiniteElement::GenericFiniteElement *elem)
          -> indexType override;

      void toParaviewAdaper(
          PointerCollection& pointers, FiniteElement::GenericFiniteElement* elem,
          vtkPlotInterface& paraviewAdapter, ParaviewSwitch control) override;
      
      auto  computeNorm(PointerCollection& pointers, FiniteElement::GenericFiniteElement* elem, FiniteElement::NormTypes type) -> prec override;


    private:

      void setTangentResidualFluss(
        PointerCollection& pointers,
        FiniteElement::GenericFiniteElement* elem,
        Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic>& stiffness,
        Eigen::Matrix<prec, Eigen::Dynamic, 1>& residual, std::vector<DegreeOfFreedom*>& Dofs);
        
      Types::VectorX<prec> getBMatrixFluss(const Types::VectorX<prec>& shapeDeriv, indexType numDofs);

      void setTangentResidualHellingerReissner(
        PointerCollection& pointers,
        FiniteElement::GenericFiniteElement* elem,
        Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic>& stiffness,
        Eigen::Matrix<prec, Eigen::Dynamic, 1>& residual, std::vector<DegreeOfFreedom*>& Dofs);

      void setTangentResidualtest(
        PointerCollection& pointers,
        FiniteElement::GenericFiniteElement* elem,
        Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic>& stiffness,
        Eigen::Matrix<prec, Eigen::Dynamic, 1>& residual, std::vector<DegreeOfFreedom*>& Dofs);


      Types::MatrixXX<prec> getBMatrixepsilon(const Types::Matrix2X<prec>& shapeDerivdisp, indexType numDofsdisp);
      Types::MatrixXX<prec> getBMatrixstress(const Types::Matrix2X<prec>& shapestress, indexType numDofsstress);
      Types::MatrixXX<prec> getBMatrixstressSkew(const Types::Matrix2X<prec>& shapestress, indexType numDofsstress);
      Types::MatrixXX<prec> getMaterial(FiniteElement::GenericFiniteElement* elem);

      
    prec E, nu;
    indexType plainstrain, meshiddisp, meshidstress, disporder, stressorder, type, mode;



  };


}
