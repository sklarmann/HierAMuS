// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include "geometry/Volumes/VolumesData.h"
#include "geometry/Volumes/VolumesRuntime.h"


namespace HierAMuS::Geometry {

VolumesRuntime::VolumesRuntime(VolumesData &data_element)
    : GeometryBaseRuntime(data_element) {}

VolumesRuntime::~VolumesRuntime() = default;

auto VolumesRuntime::getGroupType() -> const GeometryTypes & { return HierAMuS::Geometry::VolumesRuntime::type; }



const GeometryTypes VolumesRuntime::type = GeometryTypes::Volumes;

} /* namespace HierAMuS::Geometry */

