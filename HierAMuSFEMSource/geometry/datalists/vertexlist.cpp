// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause




#include <geometry/datalists/vertexlist.h>


#include <sstream>

namespace HierAMuS::Geometry {

vertexlist::vertexlist() : maxVertex(0), numberOfVertices(0) {}

vertexlist::~vertexlist() {
  this->vertIdMap.clear();
  this->vertices.clear();
  this->maxVertex = 0;
}

auto vertexlist::newElement(GeometryTypes type, indexType number) -> indexType {
  if (this->vertIdMap.find(number) == this->vertIdMap.end()) {
    indexType pos1 = this->vertices.size();

    this->vertIdMap[number] = std::make_pair(type, pos1);
    this->vertices.emplace_back();
    this->vertices[pos1].setId(number);
    if (number > this->maxVertex) {
      this->maxVertex = number;
    }
    ++this->numberOfVertices;
  }
  return number;
}

auto vertexlist::getVertex(indexType number) -> Vertex & {
  if (this->vertIdMap.find(number) == this->vertIdMap.end()) {
    std::stringstream msg;
    msg << "Requested vertex with number " << number << " does not exist!";
    throw std::runtime_error(msg.str());
  }

  return this->vertices[this->vertIdMap[number].second];
}

auto vertexlist::getItemLocalNumber(indexType number) -> Vertex * {

  return &this->vertices[number];
}


auto vertexlist::getVertexClosestTo(Types::Vector3<prec> point) -> Vertex &{
  prec minDist = std::numeric_limits<prec>::max();
  indexType minIndex = 0;
  for (auto &v : this->vertices) {
    prec dist = (v.getCoordinates() - point).norm();
    if (dist < minDist) {
      minDist = dist;
      minIndex = v.getId();
    }
  }
  return this->getVertex(minIndex);
}

} // end of namespace



