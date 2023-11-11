// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause




#include "finiteElements/ElementList.h"
#include "finiteElements/Edge.h"
#include "finiteElements/Face.h"
#include "finiteElements/Volume.h"
#include "finiteElements/FaceConstraint.h"
#include "finiteElements/VolumeConstraint.h"
#include "finiteElements/beamInterfaceElement2D.h"
#include "plot/vtkplotClassBase.h"
#include "plot/vtkplotClass.h"
#include "EquationHandler.h"
#include "geometry/GeometryData.h"
#include "solver/GenericSolutionState.h"
#include <iostream>
#include <vector>

namespace HierAMuS::FiniteElement {

ElementList::ElementList() {
}

ElementList::~ElementList() {
}

auto ElementList::requestNewElement(PointerCollection& pointers, Elementtypes type)
-> FiniteElement::GenericFiniteElement * {


  indexType number = m_elements.size();
  switch (type) {
  case Elementtypes::Edge:
    m_elements.emplace_back(std::make_shared<Edge>());
    break;
  case Elementtypes::Face:
    m_elements.emplace_back(std::make_shared<Face>());
    break;
  case Elementtypes::FaceConstraint:
    m_elements.emplace_back(std::make_shared<FaceConstraint>());
    break;
  case Elementtypes::Volume:
    m_elements.emplace_back(std::make_shared<Volume>());
    break;
  case Elementtypes::VolumeConstraint:
    m_elements.emplace_back(std::make_shared<VolumeConstraint>());
    break;
  case Elementtypes::beamInterfaceElement2D:
    m_elements.emplace_back(std::make_shared<beamInterfaceElement2D>());
    break;
    /*case Elementtypes::beamInterfaceElement3D:
    number = static_cast<indexType>(this->beamInterface3D.size());
    this->elementIndex.push_back(number);
    this->elementTypes.push_back(type);
    number = static_cast<indexType>(this->elementTypes.size()) - 1;
    this->beamInterface3D.emplace_back();
    temp = &this->beamInterface3D.back();
    temp->setId(number);
    break;
  case Elementtypes::LinearPrism:
    number = static_cast<indexType>(this->linearPrisms.size());
    this->elementIndex.push_back(number);
    this->elementTypes.push_back(type);
    number = static_cast<indexType>(this->elementTypes.size()) - 1;
    this->linearPrisms.emplace_back();
    temp = &this->linearPrisms.back();
    temp->setId(number);*/

    break;
  case Elementtypes::Generic:
    std::cout << "Error" << std::endl;
    break;
  default:
    break;
  }
  m_elements.back()->setId(number);
  return &(*m_elements.back());
}

auto ElementList::getElement(PointerCollection & pointers, indexType number)
    -> FiniteElement::GenericFiniteElement * {

  m_elements[number]->set_pointers(pointers);
  return &(*m_elements[number]);
  
}


void ElementList::setDegreesOfFreedom(PointerCollection& pointers){
  pointers.getGeometryData()->createRuntimeObjects();
  pointers.getSolutionState()->setupHistoryData(pointers);
  for (auto &i : m_elements) {
    i->set_pointers(pointers);
    i->GenericSetDegreesOfFreedom(pointers);
  }
  pointers.getEquationHandler()->update();
  pointers.getEquationHandler()->updateEquations();
  pointers.getGeometryData()->updateRuntimeObjectEquations();
  for (auto &i : m_elements) {
    i->set_pointers(pointers);
    i->GenericAdditionalOperations(pointers);
  }
  pointers.getSPDLogger().info((*pointers.getEquationHandler()));
}

auto ElementList::getNumberOfElements() -> indexType {
  return indexType(this->m_elements.size());
}

void FiniteElement::ElementList::print(PointerCollection &pointers) {
  //indexType numElems = this->getNumberOfElements();
}

} // namespace HierAMuS::FiniteElement

