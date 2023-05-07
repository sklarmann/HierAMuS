// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <geometry/Vertex.h>
#include <geometry/GeometryTypes.h>
#include <unordered_map>
#include <vector>

#include <geometry/datalists/datalistiterator.h>

namespace HierAMuS::Geometry {

class vertexlist {
  using iterator = datalistiterator<vertexlist, indexType, Vertex *>;
  using const_iterator = datalistiterator<vertexlist, indexType, const Vertex *>;

public:
  vertexlist();
  virtual ~vertexlist();

  auto newElement(GeometryTypes type, indexType number) -> indexType;
  auto getVertex(indexType number) -> Vertex &;
  auto lastElement() -> indexType { return this->maxVertex; };

  auto getNumberOfElements() -> indexType { return this->vertices.size(); };

  using vect_IT = typename std::vector<Vertex>::iterator;

  auto begin() -> iterator { return {this, 0}; };
  auto end() -> iterator { return {this, this->numberOfVertices}; };

  auto getItemLocalNumber(indexType number) -> Vertex *;

  auto getVertexClosestTo(Types::Vector3<prec> point) -> Vertex &;

private:
  indexType maxVertex;
  indexType numberOfVertices;
  std::vector<Vertex> vertices;
  std::unordered_map<indexType, std::pair<GeometryTypes, indexType>> vertIdMap;
};

} // namespace HierAMuS