// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include "geometry/GeometryData.h"
#include "geometry/GeometryBaseData.h"

#include "geometry/Volumes/LinearPrismRuntime.h"
#include <geometry/Edges/EdgesData.h>
#include <geometry/Faces/FacesData.h>
#include <geometry/VertexData.h>
#include <geometry/Volumes/LinearPrismData.h>

#include <stdexcept>
#include <types/MatrixTypes.h>

#include <vector>

#include "shapefunctions/LobattoShapes.h"

#include <iomanip>

namespace HierAMuS {
namespace Geometry {

LinearPrismData::LinearPrismData() : VolumesDataInterface() {}

LinearPrismData::~LinearPrismData() {}

auto LinearPrismData::getRuntimeObject(GeometryData &geoData)
    -> std::shared_ptr<VolumesRuntime> {
  return std::make_shared<LinearPrismRuntime>(geoData, *this);
}

const GeometryTypes &LinearPrismData::getType() { return this->type; }

void LinearPrismData::setH1Shapes(indexType meshId, indexType order,
                                  NodeTypes type) {
  FacesData *tempFace;
  for (auto i = 0; i < 5; i++) {
    tempFace = m_faces_pointers[i];
    tempFace->setH1Shapes(meshId, order, type);
  }
  if (order > 1) {
    this->setH1ShapesInternal(meshId, order, type);
  }
}

void LinearPrismData::setH1ShapesInternal(indexType meshId, indexType order,
                                          NodeTypes type) {
  if (order > 1) {
    indexType num = order - 1;
    num *= num * num;
    this->setNodeSet(meshId, num, type);
  }
}

void LinearPrismData::checkUpdateElement(EquationHandler &eqHandler,
                                         GeometryData &geoData) {
  std::cout << "Warning: checkUpdateElement not implemented in LinearPrism!"
            << std::endl;
}

const GeometryTypes LinearPrismData::type = GeometryTypes::LinearPrism;
} // namespace Geometry
} // namespace HierAMuS
