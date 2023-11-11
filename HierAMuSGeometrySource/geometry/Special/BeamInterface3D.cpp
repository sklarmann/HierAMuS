// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "GenericNodes.h"

#include <types/MatrixTypes.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <geometry/Special/BeamInterface3D.h>
#include <geometry/GeometryData.h>
#include <pointercollection/pointercollection.h>

#include <geometry/Edges/EdgesData.h>
#include <geometry/GeometryData.h>
#include <geometry/VertexData.h>
#include "geometry/Faces/FacesData.h"

#include <algorithm>
#include <iomanip>
#include <vector>

namespace HierAMuS::Geometry {

BeamInterface3D::BeamInterface3D() : A(0) {}

BeamInterface3D::~BeamInterface3D() = default;

void BeamInterface3D::setEdges(const std::vector<indexType> &edgesIn) {}

void BeamInterface3D::setFaces(std::vector<indexType> &faces) {
  this->m_faces = faces;
}

void BeamInterface3D::setBeamVertex(indexType vertex) { m_beamVertex = vertex; }

auto BeamInterface3D::getIntegrationPoints(indexType elementId) -> IntegrationPoints {
  auto intpoints = IntegrationPointsManagement::getIntegrationsPoints(elementId);
  intpoints.setType(IntegrationType::Scaled3D);
  intpoints.setNumberOfSections(this->m_faces.size());
  return intpoints;
}

void BeamInterface3D::computeWarpingShapes(PointerCollection &pointers) {
  this->createNodeShapeMapping(pointers);
  this->computeGeometry(pointers);
  auto GP = this->getIntegrationPoints(-1);
  GP.setOrder(m_warpOrder);

  for (auto i : GP) {
    auto faceIntegration = this->getFaceIntegrationPoint(i);
    auto face = pointers.getGeometryData()->getFaceData(m_faces[i.sectionNumber]);
    auto h1shapes = face->getH1Shapes(m_warpOrder, faceIntegration);
    auto Nodes = face->getH1Nodes(m_H1MeshId, m_warpOrder);

    Types::Vector3<prec> localCoor = this->getLocalCoordinate(pointers, i);
    std::cout << "Local coordinate at integrationpoint: "
              << localCoor.transpose() << std::endl;

    //Types::Vector3<prec> G1 = face->getTangent_G1(pointers, faceIntegration);
    //Types::Vector3<prec> G2 = face->getTangent_G2(pointers, faceIntegration);
    //this->A+= (G1.cross(G2)).dot(m_A1);
  }
  std::cout << "A: " << this->A << std::endl;
}

auto BeamInterface3D::getCoordinates(PointerCollection &pointers,
                                     IntegrationPoint &IntPoint)
    -> Types::Vector3<prec> {
  Types::Vector3<prec> coor;
  coor.setZero();
  auto &Face = *pointers.getGeometryData()->getFaceData(
      this->m_faces[IntPoint.sectionNumber]);
  IntegrationPoint faceInt;
  faceInt.xi = IntPoint.eta;
  faceInt.eta = IntPoint.zeta;
  coor = Face.getCoordinates(faceInt);
  prec xip = IntPoint.xi + prec(1);
  prec dL = this->m_length / prec(2) * xip;
  coor += dL * m_A1;
  return coor;
}

void BeamInterface3D::computeGeometry(PointerCollection &pointers) {
  auto face = pointers.getGeometryData()->getFaceData(m_faces[0]);
  Types::Vector3<prec> x1;
  Types::Vector3<prec> x2;
  Types::Vector3<prec> x3;
  x1 = face->getVertex(0)->getCoordinates();
  x2 = face->getVertex(1)->getCoordinates();
  x3 = face->getVertex(3)->getCoordinates();

  // compute the local basis system on the surface.
  m_A2 = x2 - x1;
  m_A3 = x3 - x1;
  m_A1 = m_A2.cross(m_A3);
  m_A1.normalize();
  m_A2 = m_A3.cross(m_A1);
  m_A2.normalize();
  m_A3 = m_A1.cross(m_A2);
  m_A3.normalize();

  // compute the projected coordinate.
  Types::Vector3<prec> beam_coordinate =
      pointers.getGeometryData()->getVertexData(m_beamVertex).getCoordinates();
  Types::Vector3<prec> temp = beam_coordinate - x1;
  m_length = temp.dot(m_A1);

  m_projectedCoordinate = beam_coordinate - m_length * m_A1;
}

void BeamInterface3D::createNodeShapeMapping(PointerCollection &pointers) {
  indexType ns = 0;
  for (auto faceId : m_faces) {
    auto face = pointers.getGeometryData()->getFaceData(faceId);
    auto Nodes = face->getH1Nodes(m_H1MeshId, m_warpOrder);
    for (auto &node : Nodes) {
      if (m_nodeShapeMapping.find(node->getId()) == m_nodeShapeMapping.end()) {
        m_nodeShapeMapping[node->getId()] = ns;
        ++ns;
      }
    }
  }
  m_numberOfWarpingShapes = ns;
}

void BeamInterface3D::setH1Shapes(PointerCollection &pointers, indexType meshId,
                                  indexType order, NodeTypes type) {
  m_H1MeshId = meshId;
  m_warpOrder = order;
  for (auto i : this->m_faces) {
    auto face = pointers.getGeometryData()->getFaceData(i);
    face->setH1Shapes(meshId, order, type);
  }
}

void BeamInterface3D::setH1ShapesBeamRot(PointerCollection &pointers,
                                         indexType meshId, NodeTypes type) {
  auto vert = pointers.getGeometryData()->getVertexData(m_beamVertex);
  vert.setNodeSet(meshId, 1, type);
}

auto BeamInterface3D::getLocalCoordinate(PointerCollection &pointers,
                                         IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  indexType currFaceNumber = m_faces[integrationPoint.sectionNumber];
  auto face = pointers.getGeometryData()->getFaceData(currFaceNumber);
  IntegrationPoint faceIntegrationPoint;
  faceIntegrationPoint.xi = integrationPoint.xi;
  faceIntegrationPoint.eta = integrationPoint.eta;
  faceIntegrationPoint.weight = integrationPoint.weight;

  Types::Vector3<prec> faceCoordinate =
      face->getCoordinates(faceIntegrationPoint);
  Types::Vector3<prec> tempCoordinate = faceCoordinate - m_projectedCoordinate;
  Types::Vector3<prec> localCoordniate;
  localCoordniate(0) = tempCoordinate.dot(m_A1);
  localCoordniate(1) = tempCoordinate.dot(m_A2);
  localCoordniate(2) = tempCoordinate.dot(m_A3);

  return localCoordniate;
}

auto BeamInterface3D::getFaceIntegrationPoint(
    IntegrationPoint &integrationPoint) -> IntegrationPoint {
  IntegrationPoint face;
  face.xi = integrationPoint.eta;
  face.eta = integrationPoint.zeta;
  return face;
}

const GeometryTypes BeamInterface3D::type = GeometryTypes::BeamInterface3D;

} // namespace HierAMuS::Geometry
