// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <materials/MaterialformulationList.h>
#include <pointercollection/pointercollection.h>


#include <materials/MaterialFormulation2.h>

#include <materials/3D/MA1_LinearElastic_Isotrop.h>
#include "materials/3D/MA2_NeoHook.h"
#include "materials/3D/MA3_SmallStrainPlasticity.h"

#include <materials/2D/MA3_2D_LinearElastic_Isotrop.h>
#include "materials/2D/MA1_2D_PlainStrain_3D.h"

#include "Homogenization/MAS1_Homogenization.h"

namespace HierAMuS::Materials {

MaterialFormulationList::MaterialFormulationList() {}

MaterialFormulationList::~MaterialFormulationList() {

}

void MaterialFormulationList::addMaterial(
    PointerCollection &pointers, indexType number,
    indexType materialFormulation, ParameterList &materialparameters) {
  while (number + 1 > this->Materials.size()) {
    this->Materials.push_back(nullptr);
  }

  switch (materialFormulation) {
  case 1:
    // Linear elastic Material
    this->Materials[number] = std::make_shared<MA1_LinearElastic_Isotrop>(&pointers);
    break;
  case 2:
    // Linear Piezo Material
    this->Materials[number] = std::make_shared<MaterialFormulation2>(&pointers);
    break;
  case 201:
    // Linear Elastic isotropic 2D material
    this->Materials[number] = std::make_shared<MA2_2D_PlainStrain_3D>(&pointers);
    break;
  case 203:
    // Linear Elastic isotropic 2D material
    this->Materials[number] = std::make_shared<MA3_2D_LinearElastic_Isotrop>(&pointers);
    break;
  case 301:
    // Linear elastic Material
    this->Materials[number] = std::make_shared<MA1_LinearElastic_Isotrop>(&pointers);
    break;
  case 302:
    // Linear elastic Material
    this->Materials[number] = std::make_shared<MA2_NeoHook>(&pointers);
    break;
  case 303:
    // Linear elastic Material
    this->Materials[number] = std::make_shared<MA3_SmallStrainPlasticity>(&pointers);
    break;
  case 400:
    this->Materials[number] = std::make_shared<MAS1_Homogenization>(&pointers);
    break;
  default:
    break;
  }
  auto temp = this->Materials[number];

  if (temp != nullptr)
    temp->readData(pointers,materialparameters);
}

std::shared_ptr<GenericMaterialFormulation>
MaterialFormulationList::getMaterial(indexType number) {
  return this->Materials[number];
}
void MaterialFormulationList::addMaterial(
    indexType number, std::shared_ptr<GenericMaterialFormulation> material) {
  while (number + 1 > this->Materials.size()) {
    this->Materials.push_back(nullptr);
  }
  this->Materials[number] = material;
}
} // namespace HierAMuS
