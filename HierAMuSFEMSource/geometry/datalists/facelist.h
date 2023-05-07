// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once


#include <geometry/LinearTriangle.h>
#include <geometry/LinearQuadrilateral.h>
#include <geometry/ScaledBoundary2D.h>
#include "geometry/QuadraticQuadrilateral.h"


#include <geometry/GeometryTypes.h>
#include <unordered_map>
#include <vector>

#include <geometry/datalists/datalistiterator.h>

namespace HierAMuS::Geometry {

class facelist {
  using iterator = datalistiterator<facelist, indexType, Faces *>;
  using const_iterator = datalistiterator<facelist, indexType, const Faces *>;

public:
  facelist();
  virtual ~facelist();

  auto newElement(GeometryTypes type, indexType number) -> indexType;
  auto getFace(indexType number) -> Faces &;
  auto lastElement() -> indexType { return this->maxFace; };

  auto getNumberOfElements() -> indexType { return this->numFaces; };

  /** @brief Returns the face number of the face which has the vertices vert1, vert2 and vert3.
   *
   *  detailed description
   *
   *  @param [in] vert1 The global number of the start vertex of the edge.
   *  @param [in] vert2 The global number of the end vertex of the edge.
   *  @param [in] vert3 The global number of the end vertex of the edge.
   *  @return The global number of the face, returns -1 if face does not exist.
   */
  auto getFaceNumberByVertexNumbers(indexType vert1, indexType vert2, indexType vert3) -> indexType;

  // Iterator functions
  auto begin() -> iterator { return {this, 0}; };
  auto end() -> iterator { return {this, this->numFaces}; };
  auto getItemLocalNumber(indexType number) -> Faces *;

private:
  indexType maxFace;
  indexType numFaces;
  std::vector<LinearTriangle> linearTriangles;
  std::vector<LinearQuadrilateral> linearQuadrilaterals;
  std::vector<QuadraticQuadrilateral> quadraticQuadrilaterals;
  std::vector<ScaledBoundary2D> scaledBoundary;

  std::unordered_map<indexType, std::pair<GeometryTypes, indexType>> faceIdMap;
};

} // namespace HierAMuS