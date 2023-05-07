// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include "materials/2D/MA1_2D_PlainStrain_3D.h"
#include "materials/GenericMaterialFormulation.h"
#include "materials/MaterialformulationList.h"
#include "pointercollection/pointercollection.h"

#include <control/ParameterList.h>


#include <string>

namespace HierAMuS::Materials {

MA2_2D_PlainStrain_3D::MA2_2D_PlainStrain_3D(PointerCollection *ptrCol)
    : GenericMaterialFormulation(ptrCol) {}

MA2_2D_PlainStrain_3D::~MA2_2D_PlainStrain_3D() = default;

void MA2_2D_PlainStrain_3D::readData(PointerCollection& pointers, ParameterList &list) {

  m_3D_matNum = list.getIndexVal("number");
}

void MA2_2D_PlainStrain_3D::getMaterialData(
  PointerCollection& pointers, MaterialTransferData &materialInOut, IntegrationPoint& ip) {

  materialInOut.materialTangent.resize(3, 3);
  materialInOut.materialTangent.setZero();
  materialInOut.stresses.resize(3);
  materialInOut.stresses.setZero();


  auto material3D =
      pointers.getMaterialFormulationList()->getMaterial(m_3D_matNum);

  MaterialTransferData material3D_in;
  material3D_in.historyData = materialInOut.historyData;
  material3D_in.strains.resize(6);
  material3D_in.strains.setZero();
  material3D_in.strains(0) = materialInOut.strains(0);
  material3D_in.strains(1) = materialInOut.strains(1);
  material3D_in.strains(3) = materialInOut.strains(2);


  material3D->getMaterialData(pointers,material3D_in,ip);

  materialInOut.materialTangent.block(0, 0, 2, 2) = material3D_in.materialTangent.block(0, 0, 2, 2);
  materialInOut.materialTangent.block(2, 0, 1, 2) = material3D_in.materialTangent.block(3, 0, 1, 2);
  materialInOut.materialTangent.block(0, 2, 2, 1) = material3D_in.materialTangent.block(0, 3, 2, 1);
  materialInOut.materialTangent(2,2) = material3D_in.materialTangent(3,3);

  materialInOut.stresses(0) = material3D_in.stresses(0);
  materialInOut.stresses(1) = material3D_in.stresses(1);
  materialInOut.stresses(2) = material3D_in.stresses(3);
}



auto MA2_2D_PlainStrain_3D::getHistoryDataStructure(PointerCollection& pointers)
-> const HistoryDataStructure & {
  auto material3D =
      pointers.getMaterialFormulationList()->getMaterial(m_3D_matNum);
  return material3D->getHistoryDataStructure(pointers);
}

auto MA2_2D_PlainStrain_3D::getInternalVariables(
  PointerCollection& pointers, MaterialTransferData &inoutData)
-> std::map<std::string, Types::VectorX<prec>> {

  auto material3D = pointers.getMaterialFormulationList()->getMaterial(m_3D_matNum);

  return material3D->getInternalVariables(pointers,inoutData);
}


} // namespace HierAMuS::Materials
