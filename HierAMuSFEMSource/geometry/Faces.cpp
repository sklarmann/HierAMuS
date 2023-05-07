// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <geometry/Faces.h>
#include <pointercollection/pointercollection.h>

#include <iomanip>

namespace HierAMuS::Geometry {

Faces::Faces() : Base() {}

Faces::~Faces() {}

const GeometryTypes &Faces::getGroupType() { return this->type; }

const GeometryTypes Faces::type = GeometryTypes::Faces;

void Faces::setAllNodeBoundaryConditionMeshId(PointerCollection &pointers,
                                         indexType meshId, indexType dof) {
  Base::setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);
}

} // namespace HierAMuS
