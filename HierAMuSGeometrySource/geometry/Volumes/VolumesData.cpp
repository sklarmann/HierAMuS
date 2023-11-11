// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include <geometry/Volumes/VolumesData.h>

#include <iomanip>

namespace HierAMuS::Geometry {

VolumesData::VolumesData() : GeometryBaseData() {}

VolumesData::~VolumesData() = default;

auto VolumesData::getGroupType() -> const GeometryTypes & { return HierAMuS::Geometry::VolumesData::type; }



const GeometryTypes VolumesData::type = GeometryTypes::Volumes;

} /* namespace HierAMuS::Geometry */

