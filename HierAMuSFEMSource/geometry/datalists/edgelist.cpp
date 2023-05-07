// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "datatypes.h"
#include <geometry/datalists/edgelist.h>


#include <sstream>


namespace HierAMuS::Geometry {

edgelist::edgelist() {
  this->maxEdges = 0;
  this->numEdges = 0;
}

edgelist::~edgelist() {
  this->linearEdges.clear();
  this->edgeIdMap.clear();
  this->maxEdges = 0;
  this->numEdges = 0;
}

auto edgelist::newElement(GeometryTypes type, indexType number) -> indexType {

  if (this->edgeIdMap.find(number) == this->edgeIdMap.end()) {
    switch (type) {
    case GeometryTypes::LinearEdge: {
      indexType pos1 = this->linearEdges.size();

      this->edgeIdMap[number] = {type, pos1};
      this->linearEdges.emplace_back();
      this->linearEdges[pos1].setId(number);

    } break;
    case GeometryTypes::QuadraticEdge: {
      indexType pos1 = this->quadraticEdges.size();

      this->edgeIdMap[number] = {type, pos1};
      this->quadraticEdges.emplace_back();
      this->quadraticEdges[pos1].setId(number);

    } break;
    default:
      std::stringstream msg;
      msg << "Trying to add an edge with the wrong type!";
      throw std::runtime_error(msg.str());
    }
    if (number > this->maxEdges) {
      this->maxEdges = number;
    }
    ++this->numEdges;
  }
  return number;
}

auto edgelist::getEdge(indexType number) -> Edges & {
  if (this->edgeIdMap.find(number) == this->edgeIdMap.end()) {
    std::stringstream msg;
    msg << "Requested Edge with number " << number << " does not exist!";
    throw std::runtime_error(msg.str());
  }
  auto TypePos = this->edgeIdMap[number];

  switch (TypePos.first) {
  case GeometryTypes::LinearEdge:
    return this->linearEdges[TypePos.second];
  case GeometryTypes::QuadraticEdge:
    return this->quadraticEdges[TypePos.second];
  default:
    std::stringstream msg;
    msg << "Requested Edge with number " << number << " does not exist!"
        << std::endl;
    throw std::runtime_error(msg.str());
  }
}

auto edgelist::getItemLocalNumber(indexType number) -> Edges * {
  indexType posA = linearEdges.size();
  indexType posB = posA + quadraticEdges.size();
  if (number < posA) {
    return &this->linearEdges[number];
  } else if (number < posB) {
    return &this->quadraticEdges[number - posA];
  }
  return nullptr;
}

auto edgelist::getEdgeNumberByVertexNumbers(indexType vert1, indexType vert2)
    -> indexType {
  for (auto &edge : this->linearEdges) {
    if (edge.hasVertices(vert1, vert2)) {
      return edge.getId();
    }
  }

  return -1;
}



} // namespace HierAMuS::Geometry
