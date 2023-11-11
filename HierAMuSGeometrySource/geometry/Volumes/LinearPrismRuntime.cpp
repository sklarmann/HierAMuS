// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "geometry/GeometryData.h"

#include "geometry/Volumes/LinearPrismData.h"
#include <geometry/Edges/EdgesData.h>
#include "geometry/Edges/EdgesRuntime.h"
#include "geometry/Edges/EdgeH1ShapesInterface.h"
#include <geometry/Faces/FacesData.h>
#include "geometry/Faces/FacesRuntime.h"
#include "geometry/Faces/FacesH1Interface.h"
#include <geometry/VertexData.h>
#include <geometry/Volumes/LinearPrismRuntime.h>

#include <stdexcept>
#include <types/MatrixTypes.h>

#include <vector>

#include "shapefunctions/LobattoShapes.h"

#include <iomanip>

namespace HierAMuS {
namespace Geometry {

LinearPrismRuntime::LinearPrismRuntime(GeometryData &geoData,
                                       LinearPrismData &data_element)
    : VolumesRuntimeDataInterface(geoData, data_element),
      m_LinearPrism_data(data_element) {}

LinearPrismRuntime::~LinearPrismRuntime() {}

const GeometryTypes &LinearPrismRuntime::getType() { return this->type; }

void LinearPrismRuntime::print(spdlog::logger &Log) {
  //
  m_LinearPrism_data.print(Log);
}

void LinearPrismRuntime::getEdgeNumbers(std::vector<indexType> &edgesOut) {
  edgesOut.assign(std::begin(this->edges), std::end(this->edges));
}

void LinearPrismRuntime::getFaceNumbers(std::vector<indexType> &facesOut) {
  facesOut.assign(std::begin(this->faces), std::end(this->faces));
}

void LinearPrismRuntime::setH1Shapes(indexType meshId, indexType order,
                                     NodeTypes type) {
  
  for (auto i = 0; i < 5; i++) {
    //auto &tempFace = *m_Faces[i];
    m_Faces[i]->getH1Face()->setH1Shapes(meshId, order, type);
  }
  if (order > 1) {
    this->setH1ShapesInternal(meshId, order, type);
  }
}

void LinearPrismRuntime::setH1ShapesInternal(indexType meshId, indexType order,
                                             NodeTypes type) {
  if (order > 1) {
    indexType num = order - 1;
    num *= num * num;
    this->setNodeSet(meshId, num, type);
  }
}

void LinearPrismRuntime::getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs,
                                   indexType meshID, indexType order) {

  for (auto i = 0; i < 6; i++) {
    auto &tempGeo = *m_Vertices[i];

    auto nodeList = tempGeo.getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }

  if (order > 1) {
    for (auto i = 0; i < 9; i++) {
      m_Edges[i]->getH1Edge()->getH1DofsInternal(Dofs, meshID, order);
    }
    for (auto i = 0; i < 5; i++) {
      m_Faces[i]->getH1Face()->getH1DofsInternal(Dofs, meshID, order);
    }
    this->getH1DofsInternal(Dofs, meshID, order);
  }
}

void LinearPrismRuntime::getH1DofsInternal(std::vector<DegreeOfFreedom *> &Dofs,
                                           indexType meshID, indexType order) {
  if (order > 1) {
    auto nodeList = this->getNodeSetNodeListMeshId(meshID);
    auto tdofs(nodeList.getDegreesOfFreedom());
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
}

void LinearPrismRuntime::getH1Shapes(indexType order,
                                     Types::VectorX<prec> &shape,
                                     Types::Matrix3X<prec> &shapeDerivative,
                                     prec xsi, prec eta, prec zeta) {

  throw std::runtime_error(
      "Error in linear prism, method getH1Shapes not implemented!");

  indexType numshapes = (order + 1);
  numshapes *= numshapes * numshapes;

  shape.resize(numshapes);
  shapeDerivative.resize(3, numshapes);

  auto [s1, ds1] = LobattoShapes::getShape(xsi, 0);
  auto [s2, ds2] = LobattoShapes::getShape(xsi, 1);

  auto [s3, ds3] = LobattoShapes::getShape(eta, 0);
  auto [s4, ds4] = LobattoShapes::getShape(eta, 1);

  auto [s5, ds5] = LobattoShapes::getShape(zeta, 0);
  auto [s6, ds6] = LobattoShapes::getShape(zeta, 1);

  shape(0) = s1 * s3 * s5;
  shapeDerivative(0, 0) = ds1 * s3 * s5;
  shapeDerivative(1, 0) = s1 * ds3 * s5;
  shapeDerivative(2, 0) = s1 * s3 * ds5;

  shape(1) = s2 * s3 * s5;
  shapeDerivative(0, 1) = ds2 * s3 * s5;
  shapeDerivative(1, 1) = s2 * ds3 * s5;
  shapeDerivative(2, 1) = s2 * s3 * ds5;

  shape(2) = s2 * s4 * s5;
  shapeDerivative(0, 2) = ds2 * s4 * s5;
  shapeDerivative(1, 2) = s2 * ds4 * s5;
  shapeDerivative(2, 2) = s2 * s4 * ds5;

  shape(3) = s1 * s4 * s5;
  shapeDerivative(0, 3) = ds1 * s4 * s5;
  shapeDerivative(1, 3) = s1 * ds4 * s5;
  shapeDerivative(2, 3) = s1 * s4 * ds5;

  shape(4) = s1 * s3 * s6;
  shapeDerivative(0, 4) = ds1 * s3 * s6;
  shapeDerivative(1, 4) = s1 * ds3 * s6;
  shapeDerivative(2, 4) = s1 * s3 * ds6;

  shape(5) = s2 * s3 * s6;
  shapeDerivative(0, 5) = ds2 * s3 * s6;
  shapeDerivative(1, 5) = s2 * ds3 * s6;
  shapeDerivative(2, 5) = s2 * s3 * ds6;

  shape(6) = s2 * s4 * s6;
  shapeDerivative(0, 6) = ds2 * s4 * s6;
  shapeDerivative(1, 6) = s2 * ds4 * s6;
  shapeDerivative(2, 6) = s2 * s4 * ds6;

  shape(7) = s1 * s4 * s6;
  shapeDerivative(0, 7) = ds1 * s4 * s6;
  shapeDerivative(1, 7) = s1 * ds4 * s6;
  shapeDerivative(2, 7) = s1 * s4 * ds6;

  if (order > 1) {

    Types::VectorX<prec> tempshape;
    Types::Matrix3X<prec> tempdshape;
    this->getH1ShapesInternal(order, tempshape, tempdshape, xsi, eta,
                              zeta);
  }
}

void LinearPrismRuntime::getH1ShapesInternal(
    indexType order, Types::VectorX<prec> &shape,
    Types::Matrix3X<prec> &shapeDerivative, prec xsi, prec eta, prec zeta) {
  if (order > 1) {
    throw std::runtime_error(
        "Higher order shapes not implemented for linear Prism element!");
  }
}

auto LinearPrismRuntime::getH1Shapes(indexType order,
                                     IntegrationPoint &IntegrationPt)
    -> H1Shapes {
  return H1Shapes(2, 3);
}

auto LinearPrismRuntime::getH1ShapesInternal(indexType order,
                                             IntegrationPoint &IntegrationPt)
    -> H1Shapes {
  return H1Shapes(2, 3);
}

auto LinearPrismRuntime::getH1Nodes(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  return std::vector<GenericNodes *>();
}

auto LinearPrismRuntime::getH1NodesInternal(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  return std::vector<GenericNodes *>();
}

auto LinearPrismRuntime::getCoordinates(IntegrationPoint &IntPoint)
    -> Types::Vector3<prec> {
  std::cout << "Warning: getCoordinates not implemented in LinearPrism!"
            << std::endl;
  return {};
}

const GeometryTypes LinearPrismRuntime::type = GeometryTypes::LinearPrism;
} // namespace Geometry
} // namespace HierAMuS
