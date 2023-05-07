// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "datatypes.h"
#include <forwarddeclaration.h>
#include <geometry/LinearEdge.h>

#include <geometry/GeometryTypes.h>
#include <types/MatrixTypes.h>

#include <array>
#include <vector>

namespace HierAMuS::Geometry {

class QuadraticEdge : public LinearEdge {
public:
  QuadraticEdge();
  ~QuadraticEdge() override;

  
  void setVerts(GeometryData &geoData,
                std::vector<indexType> &vertsIn) override;

  /**
   * @brief Changes orientation of the face (reverts the face normal vector).
   */
  void flip() override;

private:
  indexType m_centerVertex;
  static const GeometryTypes type;
};

} // namespace HierAMuS
