// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#
#include <geometry/Edges/LinearEdgeData.h>

#include <geometry/GeometryTypes.h>
#include <types/MatrixTypes.h>

#include <array>
#include <vector>

namespace HierAMuS::Geometry {

class QuadraticEdge : public LinearEdgeData {
public:
  QuadraticEdge();
  ~QuadraticEdge() override;

  
  void setVerts(GeometryData &geoData,
                std::vector<indexType> &vertsIn) override;

  /**
   * @brief Changes orientation of the face (reverts the face normal vector).
   */
  void flip() override;

  void set_geometry_pointers(GeometryData &geoData) override;

private:
  indexType m_centerVertex;
  VertexData *m_centerVertex_pointer;
  static const GeometryTypes type;
};

} // namespace HierAMuS
