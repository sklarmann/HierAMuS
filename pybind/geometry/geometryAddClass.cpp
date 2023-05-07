// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "geometryAddClass.h"

void HierAMuS::geometryAddClass::registerFunctions()
{
  this->vert.registerFunctions();
  this->geoData.registerFunctions();
  this->geoTypes.registerFunctions();
  this->geoEdges.registerFunctions();
  this->geoLinEdge.registerFunctions();
  this->base.registerFunctions();
  this->geoFaces.registerFunctions();
  this->geoVolumes.registerFunctions();
  this->geoSpecial.registerFunctions();
  this->geoScaledBoundary2D.registerFunctions();
}
