// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "ElementTypes.h"
#include "datatypes.h"

#include <iostream>
#include <vector>
#include <memory>

namespace HierAMuS {
class PointerCollection;
}

namespace HierAMuS::FiniteElement {

class GenericFiniteElement;

/**
 * @brief List containing the finite Elements.
 * Responsible for managing the finite elements in the model.
 * Contains a vector for each possible type of finite element.
 * The element data are arranged in elementIndex and elementTypes.
 * elementIndex contains the index of element n at position n in the specific
 * vector. At this position n, the vector elementTypes contains the
 * corresponding Elementtypes to identify in which vector the element is stored.
 *
 */
class ElementList {
public:
  ElementList();
  ~ElementList();

  /**
   * @brief Requests a new Element of Elementtypes type.
   * Here a new element of the specific type is generated and
   * returned as a GenericFiniteElement pointer.
   *
   * @param[in] type Type (Elementtypes) of the element which should be created.
   * @return Returns a pointer to the new finite element as a
   * GenericFiniteElement pointer.
   */
  auto requestNewElement(PointerCollection& pointers, Elementtypes type) -> GenericFiniteElement *;
  /**
   * @brief Get the finite element with the number "number".
   *
   * @param[in] number The number of the finite element which is requested.
   * @return A Pointer to the finite element with number "number" as a
   * GenericFiniteElement pointer.
   */
  auto getElement(PointerCollection & pointers, indexType number) -> GenericFiniteElement *;
  /**
   * @brief Gets the total number of finite elements.
   *
   * @return Returns the total number of finite elements as indexType.
   */
  auto getNumberOfElements() -> indexType;


  void setDegreesOfFreedom(PointerCollection& pointers);

  void print(PointerCollection &pointers);

private:

  std::vector<std::shared_ptr<GenericFiniteElement>> m_elements;
  
  //std::vector<beamInterfaceElement2D> beamInterface2D;
  //std::vector<beamInterfaceElement3D> beamInterface3D;
  //std::vector<LinearPrism> linearPrisms;
};
} // namespace HierAMuS
