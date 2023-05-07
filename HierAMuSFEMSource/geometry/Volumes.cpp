// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <geometry/Volumes.h>
#include <pointercollection/pointercollection.h>

#include <iomanip>

namespace HierAMuS::Geometry {

Volumes::Volumes() : Base() {}

Volumes::~Volumes() = default;

auto Volumes::getGroupType() -> const GeometryTypes & { return HierAMuS::Geometry::Volumes::type; }

auto Volumes::getJacobian(PointerCollection &pointers,
                          IntegrationPoint &point) -> Types::MatrixXX<prec> {
  HierAMuS::Geometry::Volumes::throwError(
      "Error when calling getJacobian for Volumes, not implemented!");
  return {};
}

const GeometryTypes Volumes::type = GeometryTypes::Volumes;

} /* namespace HierAMuS::Geometry */

