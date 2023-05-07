// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <geometry/LinearEdge.h>
#include "geometry/QuadraticEdge.h"
#include <geometry/GeometryTypes.h>
#include <unordered_map>
#include <vector>

#include <geometry/datalists/datalistiterator.h>

#include "compressedIndexToIndexList.h"

namespace HierAMuS::Geometry {

class vertexlist;

class edgelist {
  using iterator = datalistiterator<edgelist, indexType, Edges *>;
  using const_iterator = datalistiterator<edgelist, indexType, const Edges *>;

public:
  edgelist();
  virtual ~edgelist();

  auto newElement(GeometryTypes type, indexType number) -> indexType;
  auto getEdge(indexType number) -> Edges &;
  auto lastElement() -> indexType { return this->maxEdges; };

  auto getNumberOfElements() -> indexType { return this->numEdges; };

  /** @brief Returns the edge number of the edge which has the vertices vert1 and vert2.
   *
   *  detailed description
   *
   *  @param [in] vert1 The global number of the start vertex of the edge.
   *  @param [in] vert2 The global number of the end vertex of the edge.
   *  @return The global number of the edge, returns -1 if edge does not exist.
   */
  auto getEdgeNumberByVertexNumbers(indexType vert1, indexType vert2) -> indexType;

  // Iterator functions
  auto begin() -> iterator { return {this, 0}; };
  auto end() -> iterator { return {this, this->numEdges}; };
  auto getItemLocalNumber(indexType number) -> Edges *;

private:
  indexType maxEdges;
  indexType numEdges;
  std::vector<LinearEdge> linearEdges;
  std::vector<QuadraticEdge> quadraticEdges;
  std::unordered_map<indexType, std::pair<GeometryTypes, indexType>> edgeIdMap;

};

} // namespace HierAMuS