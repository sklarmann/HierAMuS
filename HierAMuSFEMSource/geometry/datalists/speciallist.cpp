// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include <geometry/datalists/speciallist.h>


#include <sstream>
#include <stdexcept>

namespace HierAMuS::Geometry {

speciallist::speciallist() {
  this->maxEdges = 0;
  this->numEdges = 0;
}

speciallist::~speciallist() {
  this->linearEdges.clear();
  this->edgeIdMap.clear();
  this->maxEdges = 0;
  this->numEdges = 0;
}

auto speciallist::newElement(GeometryTypes type, indexType number)
    -> indexType {

  if (this->edgeIdMap.find(number) == this->edgeIdMap.end()) {
    switch (type) {
    case GeometryTypes::LinearEdge: {
      indexType pos1 = this->linearEdges.size();

      this->edgeIdMap[number] = {type, pos1};
      this->linearEdges.emplace_back();
      this->linearEdges[pos1].setId(number);
      if (number > this->maxEdges) {
        this->maxEdges = number;
      }
      ++this->numEdges;

    } break;
    default:
      std::stringstream msg;
      msg << "Trying to add an edge with the wrong type!";
      throw std::runtime_error(msg.str());
    }
  }
  return number;
}

auto speciallist::getEdge(indexType number) -> Edges & {
  if (this->edgeIdMap.find(number) == this->edgeIdMap.end()) {
    std::stringstream msg;
    msg << "Requested Edge with number " << number << " does not exist!";
    throw std::runtime_error(msg.str());
  }
  indexType pos1 = this->edgeIdMap[number].second;
  GeometryTypes type = this->edgeIdMap[number].first;

  switch (type) {
  case GeometryTypes::LinearEdge:
    return this->linearEdges[pos1];
    break;
  default:
    throw std::runtime_error("Trying to get an edge with the wrong type!");
  }
}

auto speciallist::getItemLocalNumber(indexType number) -> Edges * {
  if (number < this->linearEdges.size()) {
    return &this->linearEdges[number];
  }
  return nullptr;
}

} // namespace HierAMuS
