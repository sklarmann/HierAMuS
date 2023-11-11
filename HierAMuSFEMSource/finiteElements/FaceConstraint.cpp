// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "FaceConstraint.h"
#include "geometry/GeometryBaseData.h"
#include "geometry/GeometryData.h"
#include "geometry/Faces/FacesData.h"
#include "geometry/Faces/FacesRuntime.h"
#include "geometry/Faces/FacesH1Interface.h"
#include "pointercollection/pointercollection.h"
#include "solver/GenericSolutionState.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <vector>
#include "geometry/Faces/FacesData.h"
#include "GenericNodes.h"

namespace HierAMuS::FiniteElement {
FaceConstraint::~FaceConstraint() = default;

auto FaceConstraint::getType() -> Elementtypes { return Elementtypes::Face; }
void FaceConstraint::set_pointers(PointerCollection &pointers) {}
void FaceConstraint::setVerts(std::vector<indexType> &vertsIn) {
  this->m_vertex = vertsIn[0];
}
void FaceConstraint::setFace(indexType faceIn) { this->m_face = faceIn; }
auto FaceConstraint::getVertexIds(PointerCollection &pointers)
    -> std::vector<indexType> {
  auto tempFace = pointers.getGeometryData()->getFaceData(this->m_face);
  return tempFace->getVertexNumbers();
}
auto FaceConstraint::getVertex(FaceConstraint::ptrCol &pointers,
                               indexType localNumber) -> Geometry::VertexData & {
  std::vector<indexType> vertIds;
  vertIds = this->getVertexIds(pointers);

  return pointers.getGeometryData()->getVertexData(vertIds[localNumber]);
}
auto FaceConstraint::getEdge(FaceConstraint::ptrCol &pointers,
                             indexType localNumber) -> Geometry::EdgesData & {
  auto temp = pointers.getGeometryData()->getFaceData(this->m_face);
  auto EdgeNums = temp->getEdgeNumbers();
  return pointers.getGeometryData()->getEdgeData(EdgeNums[localNumber]);
}
auto FaceConstraint::getFace(FaceConstraint::ptrCol &pointers,
                             indexType localNumber) -> Geometry::FacesData * {
  return pointers.getGeometryData()->getFaceData(this->m_face);
}

auto FaceConstraint::getLocalCoordinateSystem(PointerCollection &pointers)
    -> Types::Matrix33<prec> {

  Types::Matrix33<prec> R;
  IntegrationPoint integrationPoint;
  integrationPoint.xi = prec(0);
  integrationPoint.eta = prec(0);

  auto face = this->getFace(pointers, 0);
  Types::Vector3<prec> Normal = face->getFaceNormal();
  Types::Vector3<prec> Tangent = face->getTangent_G1(integrationPoint);
  Tangent.normalize();
  Types::Vector3<prec> A2 = Normal.cross(Tangent).normalized();
  Types::Vector3<prec> A3 = Tangent.cross(A2).normalized();

  R.col(0) = Normal;
  R.col(1) = A2;
  R.col(2) = A3;

  return R;
}

void FaceConstraint::setAllNodeBoundaryConditionMeshId(
    FaceConstraint::ptrCol &pointers, indexType meshId, indexType dof) {
  auto temp = pointers.getGeometryData()->getFaceData(this->m_face);
  temp->setAllNodeBoundaryConditionMeshId(meshId, dof);
}

auto FaceConstraint::getNumberOfVertices(PointerCollection &pointers)
    -> indexType {
  return GenericFiniteElement::getNumberOfVertices(pointers);
}
auto FaceConstraint::getNumberOfEdges(PointerCollection &pointers)
    -> indexType {
  return GenericFiniteElement::getNumberOfEdges(pointers);
}

auto FaceConstraint::getJacobian(ptrCol &pointers,
                                 IntegrationPoint &IntegrationPt)
    -> Types::MatrixXX<prec> {
  auto geomElement = pointers.getGeometryData()->getFaceData(this->m_face);
  return geomElement->getJacobian(IntegrationPt);
}


void FaceConstraint::setH1Shapes(FaceConstraint::ptrCol &pointers,
                                 indexType meshid, indexType order) {
  Geometry::FacesData *tempFace;
  tempFace = pointers.getGeometryData()->getFaceData(this->m_face);
  tempFace->setH1Shapes(meshid, order, NodeTypes::displacement);
}

void FaceConstraint::getH1Dofs(FaceConstraint::ptrCol &pointers,
                               std::vector<DegreeOfFreedom *> &Dofs,
                               indexType meshID, indexType order) {
  pointers.getGeometryData()
      ->getFaceData(this->m_face)
      ->getH1Dofs(Dofs, meshID, order);
}

auto FaceConstraint::getH1Nodes(ptrCol &pointers, indexType meshID,
                                indexType order)
    -> std::vector<GenericNodes *> {
  return pointers.getGeometryData()
      ->getFaceData(this->m_face)
      ->getH1Nodes(meshID, order);
}



auto FaceConstraint::getH1Shapes(ptrCol &pointers, indexType order,
                                 Types::MatrixXX<prec> &jacobi,
                                 IntegrationPoint &IntegrationPt)
    -> Geometry::H1Shapes {
  auto geoFace = pointers.getGeometryData()->getFaceData(this->m_face);
  auto shapes = geoFace->getH1Shapes(order, IntegrationPt);
  shapes.shapeDeriv = jacobi.inverse().transpose() * shapes.shapeDeriv;
  return shapes;
}

auto FaceConstraint::getIntegrationPoints(ptrCol &pointers)
    -> IntegrationPoints {
  auto temp = pointers.getGeometryData()->getFaceData(this->m_face);
  return temp->getIntegrationPoints(this->m_id);
}

void FaceConstraint::geometryToParaview(PointerCollection &pointers,
                                        vtkPlotInterface &paraviewAdapter,
                                        indexType mainMesh, indexType subMesh) {
  auto temp = pointers.getGeometryData()->getFaceRuntime(this->m_face);
  temp->geometryToParaview(paraviewAdapter, mainMesh, subMesh);
};

void FaceConstraint::computeWeightsParaview(PointerCollection &pointers,
                                            vtkPlotInterface &paraviewAdapter,
                                            indexType mainMesh,
                                            indexType subMesh) {
  auto temp = pointers.getGeometryData()->getFaceRuntime(this->m_face);
  temp->computeWeightsParaview(paraviewAdapter, mainMesh, subMesh);
}

void FaceConstraint::H1SolutionToParaview(PointerCollection &pointers,
                                          vtkPlotInterface &paraviewAdapter,
                                          indexType mainMesh, indexType subMesh,
                                          indexType meshId, indexType order,
                                          std::string name) {
  auto geoElem = pointers.getGeometryData()->getFaceRuntime(this->m_face);
  auto Dofs = geoElem->getH1Face()->getH1Dofs(meshId, order);
  Types::VectorX<prec> sol = pointers.getSolutionState()->getSolution(Dofs);
  geoElem->H1SolutionToParaview(paraviewAdapter, mainMesh, subMesh,
                                order, sol, name);
}

void FaceConstraint::projectDataToParaviewVertices(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {
  auto geoElem = pointers.getGeometryData()->getFaceRuntime(this->m_face);
  geoElem->projectDataToParaviewVertices(paraviewAdapter, mainMesh,
                                         subMesh, order, IntegrationPt, data,
                                         numberComponents, name);
}

auto FaceConstraint::getVertexCoordinates(ptrCol &pointers)
    -> Types::Vector3<prec> {
  return pointers.getGeometryData()->getVertexData(this->m_vertex).getCoordinates();
}

auto FaceConstraint::getFaceCoordinates(ptrCol &pointers,
                                        IntegrationPoint &IntegrationPt)
    -> Types::Vector3<prec> {
  return pointers.getGeometryData()
      ->getFaceData(this->m_face)
      ->getCoordinates(IntegrationPt);
}

void FaceConstraint::setVertexNodes(PointerCollection &pointers,
                                    indexType meshId) {
  auto &vertex = pointers.getGeometryData()->getVertexData(this->m_vertex);
  vertex.setNodeSet(meshId, 1, NodeTypes::displacement);
}

void FaceConstraint::getVertexDofs(PointerCollection &pointers,
                   std::vector<DegreeOfFreedom *> &Dofs, indexType meshID) {
  auto &vertex = pointers.getGeometryData()->getVertexData(this->m_vertex);
  std::vector<GenericNodes *> nodeVector;
  vertex.getNodes(nodeVector, meshID);
  for (auto &node : nodeVector) {
    auto addDofs = node->getDegreesOfFreedom();
    Dofs.insert(Dofs.end(), addDofs.begin(), addDofs.end());
  }
}


auto FaceConstraint::getTangentG1(PointerCollection &pointers, IntegrationPoint &IntegrationPt)
      -> Types::Vector3<prec>{
        return pointers.getGeometryData()->getFaceData(this->m_face)->getTangent_G1(IntegrationPt);
      }
  auto FaceConstraint::getTangentG2(PointerCollection &pointers, IntegrationPoint &IntegrationPt)
      -> Types::Vector3<prec>{
        
        return pointers.getGeometryData()->getFaceData(this->m_face)->getTangent_G2(IntegrationPt);
      }

} // namespace HierAMuS::FiniteElement
