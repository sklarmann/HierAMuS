// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include <geometry/Edges/EdgesRuntime.h>
#include "geometry/Edges/EdgesData.h"

#include <types/MatrixTypes.h>

#include <iomanip>

namespace HierAMuS::Geometry {


EdgesRuntime::EdgesRuntime(EdgesData &data_element) : GeometryBaseRuntime(data_element) {}

EdgesRuntime::~EdgesRuntime() = default;

auto EdgesRuntime::getGroupType() -> const GeometryTypes & { return HierAMuS::Geometry::EdgesRuntime::type; }

const GeometryTypes EdgesRuntime::type = GeometryTypes::Edges;

void EdgesRuntime::setAllNodeBoundaryConditionMeshId(indexType meshId, indexType dof) {
  GeometryBaseRuntime::setAllNodeBoundaryConditionMeshId(meshId, dof);
}

} // namespace HierAMuS
