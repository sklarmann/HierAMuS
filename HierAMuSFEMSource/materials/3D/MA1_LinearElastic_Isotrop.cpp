// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include <materials/3D/MA1_LinearElastic_Isotrop.h>

#include <control/ParameterList.h>


#include <string>

namespace HierAMuS::Materials {

MA1_LinearElastic_Isotrop::MA1_LinearElastic_Isotrop(PointerCollection *ptrCol)
    : GenericMaterialFormulation(ptrCol) {}

MA1_LinearElastic_Isotrop::~MA1_LinearElastic_Isotrop() = default;

void MA1_LinearElastic_Isotrop::readData(PointerCollection& pointers, ParameterList &list) {

  this->m_emodul = list.getPrecVal("emodul");
  this->m_nu = list.getPrecVal("nu");


  this->m_material_tangent.setZero();
  for (auto i = 0; i < 3; i++) {
    this->m_material_tangent(i, i) = prec(1) - this->m_nu;
    this->m_material_tangent(i + 3, i + 3) =
        (prec(1) - prec(2) * this->m_nu) / prec(2);
  }
  this->m_material_tangent(0, 1) = this->m_nu;
  this->m_material_tangent(1, 0) = this->m_nu;
  this->m_material_tangent(0, 2) = this->m_nu;
  this->m_material_tangent(2, 0) = this->m_nu;
  this->m_material_tangent(1, 2) = this->m_nu;
  this->m_material_tangent(2, 1) = this->m_nu;

  this->m_material_tangent *=
      this->m_emodul / (prec(1) + this->m_nu) / (prec(1) - prec(2) * this->m_nu);
}

void MA1_LinearElastic_Isotrop::getMaterialData(
  PointerCollection& pointers, MaterialTransferData &materialInOut, IntegrationPoint& ip) {

  materialInOut.stresses = this->m_material_tangent * materialInOut.strains;
  materialInOut.materialTangent = this->m_material_tangent;
}

} // namespace HierAMuS
