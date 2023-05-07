// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "ElementTypesWrapper.h"

namespace HierAMuS {
namespace FiniteElement {

void ElementtypesWrapper::registerFunctions() {
  this->temp.value("Generic", Elementtypes::Generic)
      .value("Edge", Elementtypes::Edge)
      .value("Face", Elementtypes::Face)
      .value("FaceConstraint", Elementtypes::FaceConstraint)
      .value("Volume", Elementtypes::Volume)
      .value("VolumeConstraint", Elementtypes::VolumeConstraint)
      .value("LinearPrism", Elementtypes::LinearPrism)
      .value("beamInterfaceElement2D", Elementtypes::beamInterfaceElement2D)
      .value("beamInterfaceElement3D", Elementtypes::beamInterfaceElement3D)
	;
}
} // namespace FiniteElement
}


