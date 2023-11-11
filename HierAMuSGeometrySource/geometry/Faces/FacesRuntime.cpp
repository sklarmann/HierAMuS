// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause




#include <geometry/Faces/FacesRuntime.h>
#include "geometry/Faces/FacesData.h"

#include <iomanip>

namespace HierAMuS::Geometry {

FacesRuntime::FacesRuntime(FacesData &base_element)
    : GeometryBaseRuntime(base_element) {}

FacesRuntime::~FacesRuntime() {}

const GeometryTypes &FacesRuntime::getGroupType() { return this->type; }

const GeometryTypes FacesRuntime::type = GeometryTypes::Faces;

void FacesRuntime::setAllNodeBoundaryConditionMeshId(indexType meshId, indexType dof) {
  GeometryBaseRuntime::setAllNodeBoundaryConditionMeshId(meshId, dof);
}

} // namespace HierAMuS
