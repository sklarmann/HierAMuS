// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include <geometry/Edges/EdgesData.h>

#include <types/MatrixTypes.h>

#include <iomanip>

namespace HierAMuS::Geometry {

EdgesData::EdgesData() : GeometryBaseData() {}

EdgesData::~EdgesData() = default;

auto EdgesData::getGroupType() -> const GeometryTypes & { return HierAMuS::Geometry::EdgesData::type; }

const GeometryTypes EdgesData::type = GeometryTypes::Edges;

void EdgesData::setAllNodeBoundaryConditionMeshId(indexType meshId, indexType dof) {
  GeometryBaseData::setAllNodeBoundaryConditionMeshId(meshId, dof);
}

} // namespace HierAMuS
