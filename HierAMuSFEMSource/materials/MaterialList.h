// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once
#include "datatypes.h"

#include <vector>

namespace HierAMuS {
namespace Materials {
class Material;
class MaterialList {
public:
  MaterialList();
  ~MaterialList();
  Material *getMaterial(indexType number);
  indexType getNumberOfMaterials() {
    return static_cast<indexType>(this->list.size());
  };

private:
  std::vector<Material *> list;
};

} // namespace Materials
} // namespace HierAMuS
