// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include <geometry/Faces/FacesData.h>

#include <iomanip>

namespace HierAMuS::Geometry {

FacesData::FacesData() : GeometryBaseData() {}

FacesData::~FacesData() {}

const GeometryTypes &FacesData::getGroupType() { return this->type; }

const GeometryTypes FacesData::type = GeometryTypes::Faces;

void FacesData::setAllNodeBoundaryConditionMeshId(indexType meshId, indexType dof) {
  GeometryBaseData::setAllNodeBoundaryConditionMeshId(meshId, dof);
}

} // namespace HierAMuS
