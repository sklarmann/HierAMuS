// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "geometry/LinearEdge.h"
#include "geometry/Volumes.h"
#include <geometry/datalists/volumelist.h>


#include <sstream>

namespace HierAMuS::Geometry {

volumelist::volumelist() {
  this->maxVolumes = 0;
  this->numVolumes = 0;
}

volumelist::~volumelist() {
  this->linearBricks.clear();
  this->volumeIdMap.clear();
  this->maxVolumes = 0;
  this->numVolumes = 0;
}

auto volumelist::newElement(GeometryTypes type, indexType number) -> indexType {

  if (this->volumeIdMap.find(number) == this->volumeIdMap.end()) {
    switch (type) {
    case GeometryTypes::LinearBrick: {
      indexType pos1 = this->linearBricks.size();
      this->volumeIdMap[number] = {type, pos1};
      this->linearBricks.emplace_back();
      this->linearBricks[pos1].setId(number);

    } break;
    default:
      std::stringstream msg;
      msg << "Trying to add an volume with the wrong type!";
      throw std::runtime_error(msg.str());
    }
    if (number > this->maxVolumes) {
      this->maxVolumes = number;
    }
    ++this->numVolumes;
  }
  return number;
}

auto volumelist::getVolume(indexType number) -> Volumes & {
  if (this->volumeIdMap.find(number) == this->volumeIdMap.end()) {
    std::stringstream msg;
    msg << "Requested Volume with number " << number << " does not exist!";
    throw std::runtime_error(msg.str());
  }
  auto TypePos = this->volumeIdMap[number];

  switch (TypePos.first) {
  case GeometryTypes::LinearBrick:
    return this->linearBricks[TypePos.second];
  default:
    std::stringstream msg;
    msg << "Trying to get an volume with the wrong type!";
    throw std::runtime_error(msg.str());
  }
  return this->linearBricks[0];
}

auto volumelist::getItemLocalNumber(indexType number) -> Volumes * {
  if (number < this->linearBricks.size()) {
    return &this->linearBricks[number];
  }
  return nullptr;
}

} // namespace HierAMuS::Geometry
