// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include "geometry/Base.h"
#include "geometry/GeometryTypes.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"


#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <equations/NodeSet.h>
#include <geometry/GeometryData.h>
#include <geometry/QuadraticEdge.h>
#include <geometry/Vertex.h>
#include <pointercollection/pointercollection.h>
#include "geometry/GeometryData.h"

#include <equations/DegreeOfFreedom.h>

#include <loads/LoadList.h>

#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

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
    LinearEdge::setVerts(geoData, linearVerts);
    m_centerVertex = vertsIn[2];

  }
}

void QuadraticEdge::flip()
{
  LinearEdge::flip();
}


} // namespace HierAMuS::Geometry
