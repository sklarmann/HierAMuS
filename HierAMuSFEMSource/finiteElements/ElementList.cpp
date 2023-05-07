// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include <forwarddeclaration.h>

#include <finiteElements/ElementList.h>
#include <finiteElements/GenericFiniteElement.h>
#include "equations/EquationHandler.h"
#include <geometry/GeometryData.h>
#include <iostream>
#include <vector>

namespace HierAMuS::FiniteElement {

ElementList::ElementList() {
}

ElementList::~ElementList() {
}

auto ElementList::requestNewElement(PointerCollection& pointers, Elementtypes type)
-> FiniteElement::GenericFiniteElement * {

  FiniteElement::GenericFiniteElement *temp = nullptr;

  indexType number;
  switch (type) {
  case Elementtypes::Edge:
    number = static_cast<indexType>(this->edgeElements.size());
    this->elementIndex.push_back(number);
    this->elementTypes.push_back(type);
    number = static_cast<indexType>(this->elementTypes.size()) - 1;
    this->edgeElements.emplace_back();
    temp = &this->edgeElements.back();
    temp->setId(number);
    break;
  case Elementtypes::Face:
    number = static_cast<indexType>(this->faceElements.size());
    this->elementIndex.push_back(number);
    this->elementTypes.push_back(type);
    number = static_cast<indexType>(this->elementTypes.size()) - 1;
    this->faceElements.emplace_back();
    temp = &this->faceElements.back();
    temp->setId(number);
    break;
  case Elementtypes::FaceConstraint:
    number = static_cast<indexType>(this->faceConstraintElements.size());
    this->elementIndex.push_back(number);
    this->elementTypes.push_back(type);
    number = static_cast<indexType>(this->elementTypes.size()) - 1;
    this->faceConstraintElements.emplace_back();
    temp = &this->faceConstraintElements.back();
    temp->setId(number);
    break;
  case Elementtypes::Volume:
    number = static_cast<indexType>(this->volumeElements.size());
    this->elementIndex.push_back(number);
    this->elementTypes.push_back(type);
    number = static_cast<indexType>(this->elementTypes.size()) - 1;
    this->volumeElements.emplace_back();
    temp = &this->volumeElements.back();
    temp->setId(number);
    break;
  case Elementtypes::VolumeConstraint:
    number = static_cast<indexType>(this->volumeConstraintElements.size());
    this->elementIndex.push_back(number);
    this->elementTypes.push_back(type);
    number = static_cast<indexType>(this->elementTypes.size()) - 1;
    this->volumeConstraintElements.emplace_back();
    temp = &this->volumeConstraintElements.back();
    temp->setId(number);
    break;
  case Elementtypes::beamInterfaceElement2D:
    number = static_cast<indexType>(this->beamInterface2D.size());
    this->elementIndex.push_back(number);
    this->elementTypes.push_back(type);
    number = static_cast<indexType>(this->elementTypes.size()) - 1;
    this->beamInterface2D.emplace_back();
    temp = &this->beamInterface2D.back();
    temp->setId(number);
    break;
  case Elementtypes::beamInterfaceElement3D:
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
    temp->setId(number);

    break;
  case Elementtypes::Generic:
    std::cout << "Error" << std::endl;
    break;
  default:
    break;
  }
  return temp;
}

auto ElementList::getElement(indexType number)
    -> FiniteElement::GenericFiniteElement * {
  if (static_cast<std::size_t>(number) >= this->elementIndex.size()) {
    // TODO throw exception
  }

  switch (this->elementTypes[number]) {
  case Elementtypes::Edge:
    return &this->edgeElements[this->elementIndex[number]];
    break;
  case Elementtypes::Face:
    return &this->faceElements[this->elementIndex[number]];
    break;
  case Elementtypes::FaceConstraint:
    return &this->faceConstraintElements[this->elementIndex[number]];
    break;
  case Elementtypes::Volume:
    return &this->volumeElements[this->elementIndex[number]];
    break;
  case Elementtypes::VolumeConstraint:
    return &this->volumeConstraintElements[this->elementIndex[number]];
    break;
  case Elementtypes::beamInterfaceElement2D:
    return &this->beamInterface2D[this->elementIndex[number]];
    break;
  case Elementtypes::beamInterfaceElement3D:
    return &this->beamInterface3D[this->elementIndex[number]];
    break;
  case Elementtypes::LinearPrism:
    return &this->linearPrisms[this->elementIndex[number]];
    break;
  default:
    throw std::runtime_error(
        "Error in elementlist getElement. Requested element not in list!");
    return nullptr;
  }
  return nullptr;
}


void ElementList::setDegreesOfFreedom(PointerCollection& pointers){
  indexType numElems = this->getNumberOfElements();
  pointers.getSolutionState()->setupHistoryData(pointers);
  for(auto i=0; i<numElems; i++){
    this->getElement(i)->GenericSetDegreesOfFreedom(pointers);
  }
  pointers.getEquationHandler()->update();
  pointers.getEquationHandler()->updateEquations();
  for(auto i=0; i<numElems; i++){
    this->getElement(i)->GenericAdditionalOperations(pointers);
  }
}

} // namespace HierAMuS::FiniteElement

