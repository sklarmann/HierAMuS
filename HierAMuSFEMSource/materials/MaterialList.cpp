// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include <materials/MaterialList.h>

#include <materials/Material.h>

namespace HierAMuS::Materials {

MaterialList::MaterialList() {  }

MaterialList::~MaterialList() {
  typename std::vector<Material *>::iterator it;
  it = this->list.begin();

  while (it != this->list.end()) {
    delete *it;
    ++it;
  }
}

Material *MaterialList::getMaterial(indexType number) {
  auto currSize = static_cast<indexType>(this->list.size());
  if (number < currSize) {
    return this->list[number];
  } else {
    for (auto i = currSize; i < number + 1; ++i) {
      this->list.emplace_back();
      this->list.back() = new Material(i);
    }
    return this->list[number];
  }
}

} // namespace HierAMuS
