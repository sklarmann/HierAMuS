// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include "geometry/GeometryBaseData.h"
#include "geometry/GeometryTypes.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

#include <geometry/GeometryData.h>
#include <geometry/Edges/QuadraticEdge.h>
#include <geometry/VertexData.h>
#include "geometry/GeometryData.h"


#include "LoadList.h"

#include <shapefunctions/LagrangeShape.h>
#include <stdexcept>
#include <types/MatrixTypes.h>

#include <iomanip>
#include <sstream>

namespace HierAMuS::Geometry {
using std::vector;

QuadraticEdge::QuadraticEdge() {

}

QuadraticEdge::~QuadraticEdge() = default;

void QuadraticEdge::setVerts(GeometryData &geoData,
                             std::vector<indexType> &vertsIn)
{
  if (vertsIn.size() == 3) {
    std::vector<indexType> linearVerts = {vertsIn[0], vertsIn[1]};
    LinearEdgeData::setVerts(geoData, linearVerts);
    m_centerVertex = vertsIn[2];

  }
}

void QuadraticEdge::flip()
{
  LinearEdgeData::flip();
}

void QuadraticEdge::set_geometry_pointers(GeometryData &geoData) {
  LinearEdgeData::set_geometry_pointers(geoData);
  m_centerVertex_pointer = &geoData.getVertexData(m_centerVertex);
}

} // namespace HierAMuS::Geometry
