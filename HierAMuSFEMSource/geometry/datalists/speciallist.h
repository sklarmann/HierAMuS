// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <geometry/LinearEdge.h>
#include <geometry/GeometryTypes.h>
#include <unordered_map>
#include <vector>

#include <geometry/datalists/datalistiterator.h>

namespace HierAMuS::Geometry {

class speciallist {
  using iterator = datalistiterator<speciallist, indexType, Edges *>;
  using const_iterator = datalistiterator<speciallist, indexType, const Edges *>;

public:
  speciallist();
  virtual ~speciallist();

  auto newElement(GeometryTypes type, indexType number) -> indexType;
  auto getEdge(indexType number) -> Edges &;
  auto lastElement() -> indexType { return this->maxEdges; };

  auto getNumberOfElements() -> indexType { return this->numEdges; };

  // Iterator functions
  auto begin() -> iterator { return {this, 0}; };
  auto end() -> iterator { return {this, this->numEdges}; };
  auto getItemLocalNumber(indexType number) -> Edges *;

private:
  indexType maxEdges;
  indexType numEdges;
  std::vector<LinearEdge> linearEdges;
  std::unordered_map<indexType, std::pair<GeometryTypes, indexType>> edgeIdMap;
};

} // namespace HierAMuS