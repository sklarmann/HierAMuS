// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once

#include <forwarddeclaration.h>

#include <elementFormulations/GenericElementFormulation.h>
#include <Eigen/Dense>

namespace HierAMuS::Elementformulations {


  template<typename prec, typename indexType>
  class Hu_Washizu3DBeam : public GenericElementFormulation<prec,indexType> {
    public:
       Hu_Washizu3DBeam(PointerCollection<prec, indexType> *ptrCol);
      ~Hu_Washizu3DBeam();
      void readData(stringCommandHandler &Command, Userconstants<prec> *ucons);
      void setDegreesOfFreedom(FiniteElement::GenericFiniteElement<prec, indexType> *elem);
      void setTangentResidual(FiniteElement::GenericFiniteElement<prec, indexType> *elem
          , Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness
          , Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual
          , std::vector<DegreeOfFreedom<prec, indexType>*> &Dofs);
      
      // plot
      void projectedToParaview(PointerCollection<prec, indexType>& pointers, FiniteElement::GenericFiniteElement<prec, indexType> *elem, managementClass &paraviewAdapter);
    private:
      indexType meshIdDisp, meshIdPiezo, intOrderDisp, intOrderPiezo;

  };

}
