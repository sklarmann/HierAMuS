// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include <materials/2D/MA3_2D_LinearElastic_Isotrop.h>

#include <control/ParameterList.h>


#include <string>

namespace HierAMuS::Materials {

MA3_2D_LinearElastic_Isotrop::MA3_2D_LinearElastic_Isotrop(PointerCollection *ptrCol)
    : GenericMaterialFormulation(ptrCol) {}

MA3_2D_LinearElastic_Isotrop::~MA3_2D_LinearElastic_Isotrop() = default;

void MA3_2D_LinearElastic_Isotrop::readData(PointerCollection& pointers, ParameterList &list) {

  this->m_emodul = list.getPrecVal("emodul");
  this->m_nu = list.getPrecVal("nu");
  this->m_thick = list.getPrecVal("thickness");

  indexType pp = list.getIndexVal("plainstrain");

  pp == 1 ? this->m_plain_strain = true : this->m_plain_strain = false;

  if (this->m_plain_strain) {
    this->m_material_tangent.setZero();
    prec fac;
    fac = this->m_emodul / (prec(1) + this->m_nu) / (prec(1) - prec(2) * this->m_nu);
    this->m_material_tangent(0, 0) = (prec) 1 - this->m_nu;
    this->m_material_tangent(0, 1) = this->m_nu;
    this->m_material_tangent(1, 1) = (prec) 1 - this->m_nu;
    this->m_material_tangent(1, 0) = this->m_nu;
    this->m_material_tangent *= fac;
    this->m_material_tangent(2, 2) = this->m_emodul / (prec(1) + this->m_nu) / prec(2);

  } else {
    this->m_material_tangent.setZero();
    prec fac = this->m_emodul / (prec(1) - this->m_nu * this->m_nu);
    this->m_material_tangent(0, 0) = prec(1);
    this->m_material_tangent(0, 1) = this->m_nu;
    this->m_material_tangent(1, 1) = prec(1);
    this->m_material_tangent(1, 0) = this->m_nu;
    this->m_material_tangent *= fac;
    this->m_material_tangent(2, 2) = this->m_emodul / (prec(1) + this->m_nu) / prec(2);
    this->m_material_tangent *= this->m_thick;
  }
}

void MA3_2D_LinearElastic_Isotrop::getMaterialData(
  PointerCollection& pointers, MaterialTransferData &materialInOut, IntegrationPoint& ip) {

  materialInOut.stresses = this->m_material_tangent * materialInOut.strains;
  materialInOut.materialTangent = this->m_material_tangent;
}

} // namespace HierAMuS
