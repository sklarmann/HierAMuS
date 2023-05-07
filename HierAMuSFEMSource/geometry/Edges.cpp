// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <geometry/Edges.h>
#include <pointercollection/pointercollection.h>

#include <types/MatrixTypes.h>

#include <iomanip>

namespace HierAMuS::Geometry {

Edges::Edges() : Base() {}

Edges::~Edges() = default;

auto Edges::getGroupType() -> const GeometryTypes & { return HierAMuS::Geometry::Edges::type; }

const GeometryTypes Edges::type = GeometryTypes::Edges;

void Edges::setAllNodeBoundaryConditionMeshId(PointerCollection &pointers,
                                              indexType meshId, indexType dof) {
  Base::setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);
}

} // namespace HierAMuS
