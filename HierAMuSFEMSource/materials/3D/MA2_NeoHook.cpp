// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "MatrixTypes.h"

#include <materials/3D/MA2_NeoHook.h>

#include <control/ParameterList.h>


#include <string>

namespace HierAMuS::Materials {

MA2_NeoHook::MA2_NeoHook(PointerCollection *ptrCol)
    : GenericMaterialFormulation(ptrCol) {}

MA2_NeoHook::~MA2_NeoHook() = default;

void MA2_NeoHook::readData(PointerCollection& pointers, ParameterList &list) {
  m_G = list.getPrecVal("G");
  m_Lambda = list.getPrecVal("Lambda");

}

void MA2_NeoHook::getMaterialData(
  PointerCollection& pointers, MaterialTransferData &materialInOut, IntegrationPoint& ip) {

  materialInOut.materialTangent.resize(6, 6);
  materialInOut.materialTangent.setZero();
  materialInOut.stresses.resize(6);
  materialInOut.stresses.setZero();

  Types::Matrix33<prec> C;
  C=Types::Matrix33<prec>::Identity();
  for(auto i=0;i<3;i++){
    C(i,i) += materialInOut.strains(i)*prec(2);
  }
  C(1,0) = materialInOut.strains(3);
  C(2,0) = materialInOut.strains(4);
  C(2,1) = materialInOut.strains(5);

  C(0,1) = C(1,0);
  C(0,2) = C(2,0);
  C(1,2) = C(2,1);

  prec detC = C.determinant();

  Types::Matrix33<prec> Cinv = C.inverse();
  Types::Vector6<prec> CI;
  CI(0) = Cinv(0,0);
  CI(1) = Cinv(1,1);
  CI(2) = Cinv(2,2);
  CI(3) = Cinv(0,1);
  CI(4) = Cinv(0,2);
  CI(5) = Cinv(1,2);

  prec fact1 = prec(0.5) * m_Lambda*(detC-prec(1))-m_G;
  Types::Matrix33<prec> sig;
  for(auto i=0;i<3;i++){
    for(auto j=0;j<3;j++){
      sig(i,j) =  fact1*Cinv(i,j);
    }
    sig(i,i) += m_G;
  }

  materialInOut.stresses(0) = sig(0,0);
  materialInOut.stresses(1) = sig(1,1);
  materialInOut.stresses(2) = sig(2,2);
  materialInOut.stresses(3) = sig(0,1);
  materialInOut.stresses(4) = sig(0,2);
  materialInOut.stresses(5) = sig(1,2);


  materialInOut.materialTangent(0,0) = CI(0)*CI(0);
  materialInOut.materialTangent(0,1) = CI(3)*CI(3);
  materialInOut.materialTangent(0,2) = CI(4)*CI(4);
  materialInOut.materialTangent(0,3) = CI(0)*CI(3);
  materialInOut.materialTangent(0,4) = CI(0)*CI(4);
  materialInOut.materialTangent(0,5) = CI(3)*CI(4);

  materialInOut.materialTangent(1,1) = CI(1)*CI(1);
  materialInOut.materialTangent(1,2) = CI(1)*CI(5);
  materialInOut.materialTangent(1,3) = CI(1)*CI(3);
  materialInOut.materialTangent(1,4) = CI(3)*CI(5);
  materialInOut.materialTangent(1,5) = CI(1)*CI(5);

  materialInOut.materialTangent(2,2) = CI(2)*CI(2);
  materialInOut.materialTangent(2,3) = CI(4)*CI(5);
  materialInOut.materialTangent(2,4) = CI(2)*CI(4);
  materialInOut.materialTangent(2,5) = CI(2)*CI(5);

  materialInOut.materialTangent(3,3) = prec(0.5)*(CI(0)*CI(1) + CI(3)*CI(3));
  materialInOut.materialTangent(3,4) = prec(0.5)*(CI(0)*CI(5) + CI(3)*CI(4));
  materialInOut.materialTangent(3,5) = prec(0.5)*(CI(3)*CI(5) + CI(1)*CI(4));

  materialInOut.materialTangent(4,4) = prec(0.5)*(CI(0)*CI(2) + CI(4)*CI(4));
  materialInOut.materialTangent(4,5) = prec(0.5)*(CI(2)*CI(3) + CI(4)*CI(5));

  materialInOut.materialTangent(5,5) = prec(0.5)*(CI(1)*CI(2) + CI(5)*CI(5));

  fact1 *= -prec(2);
  prec fact2 = detC*m_Lambda;
  for(auto i=0;i<6;i++){
    for(auto j=i;j<6;j++){
      materialInOut.materialTangent(i,j) = fact1*materialInOut.materialTangent(i,j) + fact2*CI(i)*CI(j);
      materialInOut.materialTangent(j,i) = materialInOut.materialTangent(i,j);
    }
  }
}

} // namespace HierAMuS
