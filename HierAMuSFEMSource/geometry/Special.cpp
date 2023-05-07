// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <geometry/Special.h>
#include <pointercollection/pointercollection.h>

#include <iomanip>

namespace HierAMuS::Geometry {

Special::Special() : Base() {}

Special::~Special() = default;

auto Special::getGroupType() -> const GeometryTypes & {
  return HierAMuS::Geometry::Special::type;
  ;
}

const GeometryTypes Special::type = GeometryTypes::Special;

} /* namespace HierAMuS::Geometry */

