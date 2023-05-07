// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "forwarddeclaration.h"
#include "geometry/Volumes.h"
#include <geometry/LinearBrick.h>
#include <geometry/GeometryTypes.h>
#include <unordered_map>
#include <vector>

#include <geometry/datalists/datalistiterator.h>

namespace HierAMuS::Geometry {

class volumelist {
  using iterator = datalistiterator<volumelist, indexType, Volumes *>;
  using const_iterator = datalistiterator<volumelist, indexType, const Volumes *>;

public:
  volumelist();
  virtual ~volumelist();

  auto newElement(GeometryTypes type, indexType number) -> indexType;
  auto getVolume(indexType number) -> Volumes &;
  auto lastElement() -> indexType { return this->maxVolumes; };

  auto getNumberOfElements() -> indexType { return this->numVolumes; };

  // Iterator functions
  auto begin() -> iterator { return {this, 0}; };
  auto end() -> iterator { return {this, this->numVolumes}; };
  auto getItemLocalNumber(indexType number) -> Volumes *;

private:
  indexType maxVolumes;
  indexType numVolumes;
  std::vector<LinearBrick> linearBricks;
  std::unordered_map<indexType, std::pair<GeometryTypes, indexType>> volumeIdMap;
};

} // namespace HierAMuS