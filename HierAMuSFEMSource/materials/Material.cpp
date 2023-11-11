// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause




#include "materials/Material.h"

#include "materials/MaterialformulationList.h"
#include "materials/ElementformulationList.h"

#include "materials/GenericMaterialFormulation.h"
#include "elementFormulations/GenericElementFormulation.h"



#include "pointercollection/pointercollection.h"


#include <stdexcept>

namespace HierAMuS::Materials {

Material::Material(indexType matNumIn) : id(matNumIn) {
}

Material::~Material() { }

void Material::setElementForumaltion(indexType elementNumber) {
  this->elnum = elementNumber;
}

void Material::setMaterialFormulation(indexType materialNumber) {

  this->matNum = materialNumber;
}

auto Material::getMaterialFormulation(PointerCollection& pointers) -> std::shared_ptr<GenericMaterialFormulation>
{

  return pointers.getMaterialFormulationList()->getMaterial(this->matNum);
}

auto Material::getElementFormulation(
  PointerCollection& pointers) -> std::shared_ptr<Elementformulations::GenericElementFormulation>
{
  return pointers.getElementFormulationList()->getElementFormulation(
      this->elnum);
}

} /* namespace HierAMuS */
