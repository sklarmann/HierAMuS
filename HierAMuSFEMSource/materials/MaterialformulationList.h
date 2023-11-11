// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "datatypes.h"

#include <vector>
#include <memory>

namespace HierAMuS {
class PointerCollection;
class ParameterList;
namespace Materials {
class GenericMaterialFormulation;
class MaterialFormulationList {
public:
  MaterialFormulationList();
  ~MaterialFormulationList();
  void addMaterial(PointerCollection &pointers, indexType number,
                   indexType materialFormulation,
                   ParameterList &materialparameters);
  void addMaterial(indexType number,
                   std::shared_ptr<GenericMaterialFormulation> material);
  std::shared_ptr<GenericMaterialFormulation> getMaterial(indexType number);
  indexType getLastMaterialNumber() { return this->Materials.size() - 1; };

private:
  std::vector<std::shared_ptr<GenericMaterialFormulation>> Materials;
};
} // namespace Materials
} // namespace HierAMuS