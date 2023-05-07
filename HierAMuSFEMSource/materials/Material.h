// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once
#include <datatypes.h>

#include <forwarddeclaration.h>
#include <memory>

namespace HierAMuS::Materials {
class Material {
public:
  Material(indexType matNumIn);
  ~Material();
  void setElementForumaltion(indexType elementNumber);
  void setMaterialFormulation(indexType materialNumber);
  auto getElementFormulation(
    PointerCollection& pointers) -> std::shared_ptr<Elementformulations::GenericElementFormulation>;
  auto getMaterialFormulation(PointerCollection& pointers) -> std::shared_ptr<GenericMaterialFormulation>;
  auto getNumber() -> indexType { return this->id; };

private:
  indexType matNum, elnum, id;
};

} /* namespace HierAMuS */
