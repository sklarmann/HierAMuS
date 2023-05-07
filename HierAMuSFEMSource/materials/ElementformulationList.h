// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <forwarddeclaration.h>

#include <vector>
#include <elementFormulations/GenericElementFormulation.h>
#include <control/ParameterList.h>
#include <memory>

namespace HierAMuS::Materials {

class ElementFormulationList {
public:
  ElementFormulationList();
  ~ElementFormulationList();
  void addElementFormulation(PointerCollection &pointers, indexType number,
                             indexType elementFormulation,
                             ParameterList &elementparameters);
  void addElementFormulation(indexType number, std::shared_ptr<Elementformulations::GenericElementFormulation> element);
  std::shared_ptr<Elementformulations::GenericElementFormulation> 
  getElementFormulation(indexType number);
  indexType getLastElementFormulationNumber() {
    return this->Elements.size() - 1;
  };

private:
  std::vector<std::shared_ptr<Elementformulations::GenericElementFormulation>> Elements;
};
} // namespace HierAMuS