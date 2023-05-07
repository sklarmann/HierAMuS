// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "FaceConstraint.h"
#include "equations/Nodetypes.h"
#include "geometry/Base.h"
#include "pointercollection/pointercollection.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <vector>

namespace HierAMuS::FiniteElement {
FaceConstraint::~FaceConstraint() = default;

auto FaceConstraint::getType() -> Elementtypes { return Elementtypes::Face; }
void FaceConstraint::setVerts(std::vector<indexType> &vertsIn) {
  this->vertex = vertsIn[0];
}
void FaceConstraint::setFace(indexType faceIn) { Face::setFace(faceIn); }
auto FaceConstraint::getVertexIds(PointerCollection &pointers)
    -> std::vector<indexType> {
  auto tempFace = pointers.getGeometryData()->getFace(this->face);
  std::vector<indexType> vertIds;
  tempFace->getVerts(vertIds);
  return vertIds;
}
auto FaceConstraint::getVertex(FaceConstraint::ptrCol &pointers,
                               indexType localNumber) -> Geometry::Vertex & {
  std::vector<indexType> vertIds;
  vertIds = this->getVertexIds(pointers);

  return pointers.getGeometryData()->getVertex(vertIds[localNumber]);
}
auto FaceConstraint::getEdge(FaceConstraint::ptrCol &pointers,
                             indexType localNumber) -> Geometry::Edges & {
  Geometry::Base *temp;
  temp = pointers.getGeometryData()->getFace(this->face);
  std::vector<indexType> EdgeNums;
  temp->getEdges(EdgeNums);
  return pointers.getGeometryData()->getEdge(EdgeNums[localNumber]);
}
auto FaceConstraint::getFace(FaceConstraint::ptrCol &pointers,
                             indexType localNumber) -> Geometry::Faces * {
  return pointers.getGeometryData()->getFace(this->face);
}

auto FaceConstraint::getLocalCoordinateSystem(PointerCollection &pointers)
    -> Types::Matrix33<prec> {

  Types::Matrix33<prec> R;
  IntegrationPoint integrationPoint;
  integrationPoint.xi = prec(0);
  integrationPoint.eta = prec(0);

  auto face = this->getFace(pointers, 0);
  Types::Vector3<prec> Normal = face->getFaceNormal(pointers);
  Types::Vector3<prec> Tangent = face->getTangent_G1(pointers, integrationPoint);
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
  auto temp = pointers.getGeometryData()->getFace(this->face);
  temp->setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);
  // GenericFiniteElement::setAllNodeBoundaryConditionMeshId(pointers, meshId,
  // dof);
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
  auto geomElement = pointers.getGeometryData()->getFace(this->face);
  return geomElement->getJacobian(pointers, IntegrationPt);
}
void FaceConstraint::getJacobian(FaceConstraint::ptrCol &pointers,
                                 Types::Matrix22<prec> &jacobi, prec xsi,
                                 prec eta) {
  GenericFiniteElement::getJacobian(pointers, jacobi, xsi, eta);
}

void FaceConstraint::setH1Shapes(FaceConstraint::ptrCol &pointers,
                                 indexType meshid, indexType order) {
  Geometry::Faces *tempFace;
  tempFace = pointers.getGeometryData()->getFace(this->face);
  tempFace->setH1Shapes(pointers, meshid, order, NodeTypes::displacement);
}

void FaceConstraint::getH1Dofs(FaceConstraint::ptrCol &pointers,
                               std::vector<DegreeOfFreedom *> &Dofs,
                               indexType meshID, indexType order) {
  pointers.getGeometryData()
      ->getFace(this->face)
      ->getH1Dofs(pointers, Dofs, meshID, order);
}

auto FaceConstraint::getH1Nodes(ptrCol &pointers, indexType meshID,
                                indexType order)
    -> std::vector<GenericNodes *> {
  return pointers.getGeometryData()
      ->getFace(this->face)
      ->getH1Nodes(pointers, meshID, order);
}

void FaceConstraint::getH1Shapes(ptrCol &pointers, indexType order,
                                 Types::MatrixXX<prec> &jacobi,
                                 Types::VectorX<prec> &shape,
                                 Types::MatrixXX<prec> &shapeDerivative,
                                 IntegrationPoint &IntegrationPt) {
  auto geomElement = pointers.getGeometryData()->getFace(this->face);
  auto shapes = geomElement->getH1Shapes(pointers, order, IntegrationPt);
  shape = shapes.shapes;
  shapeDerivative = shapes.shapeDeriv;
  shapeDerivative = jacobi.inverse().transpose() * shapeDerivative;
}
void FaceConstraint::getH1Shapes(FaceConstraint::ptrCol &pointers,
                                 indexType order, Types::Matrix22<prec> &jacobi,
                                 Types::VectorX<prec> &shape,
                                 Types::Matrix2X<prec> &dshape, prec xi,
                                 prec eta) {
  auto geoFace = pointers.getGeometryData()->getFace(this->face);
  geoFace->getH1Shapes(pointers, order, shape, dshape, xi, eta);
  dshape = jacobi.inverse().transpose() * dshape;
}

auto FaceConstraint::getH1Shapes(ptrCol &pointers, indexType order,
                                 Types::MatrixXX<prec> &jacobi,
                                 IntegrationPoint &IntegrationPt)
    -> Geometry::H1Shapes {
  auto geoFace = pointers.getGeometryData()->getFace(this->face);
  auto shapes = geoFace->getH1Shapes(pointers, order, IntegrationPt);
  shapes.shapeDeriv = jacobi.inverse().transpose() * shapes.shapeDeriv;
  return shapes;
}

auto FaceConstraint::getIntegrationPoints(ptrCol &pointers)
    -> IntegrationPoints {
  return Face::getIntegrationPoints(pointers);
}

void FaceConstraint::geometryToParaview(PointerCollection &pointers,
                                        vtkPlotInterface &paraviewAdapter,
                                        indexType mainMesh, indexType subMesh) {
  auto temp = pointers.getGeometryData()->getFace(this->face);
  temp->geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
};

void FaceConstraint::computeWeightsParaview(PointerCollection &pointers,
                                            vtkPlotInterface &paraviewAdapter,
                                            indexType mainMesh,
                                            indexType subMesh) {
  auto temp = pointers.getGeometryData()->getFace(this->face);
  temp->computeWeightsParaview(pointers, paraviewAdapter, mainMesh, subMesh);
}

void FaceConstraint::H1SolutionToParaview(PointerCollection &pointers,
                                          vtkPlotInterface &paraviewAdapter,
                                          indexType mainMesh, indexType subMesh,
                                          indexType meshId, indexType order,
                                          std::string name) {
  auto geoElem = pointers.getGeometryData()->getFace(this->face);
  geoElem->H1SolutionToParaview(pointers, paraviewAdapter, mainMesh, subMesh,
                                meshId, order, name);
}

void FaceConstraint::projectDataToParaviewVertices(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {
  auto geoElem = pointers.getGeometryData()->getFace(this->face);
  geoElem->projectDataToParaviewVertices(pointers, paraviewAdapter, mainMesh,
                                         subMesh, order, IntegrationPt, data,
                                         numberComponents, name);
}

auto FaceConstraint::getVertexCoordinates(ptrCol &pointers)
    -> Types::Vector3<prec> {
  return pointers.getGeometryData()->getVertex(this->vertex).getCoordinates();
}

auto FaceConstraint::getFaceCoordinates(ptrCol &pointers,
                                        IntegrationPoint &IntegrationPt)
    -> Types::Vector3<prec> {
  return pointers.getGeometryData()
      ->getFace(this->face)
      ->getCoordinates(pointers, IntegrationPt);
}

void FaceConstraint::setVertexNodes(PointerCollection &pointers,
                                    indexType meshId) {
  auto &vertex = pointers.getGeometryData()->getVertex(this->vertex);
  vertex.setNodeSet(pointers, meshId, 1, NodeTypes::displacement);
}

void FaceConstraint::getVertexDofs(PointerCollection &pointers,
                   std::vector<DegreeOfFreedom *> &Dofs, indexType meshID) {
  auto &vertex = pointers.getGeometryData()->getVertex(this->vertex);
  std::vector<GenericNodes *> nodeVector;
  vertex.getNodes(pointers, nodeVector, meshID);
  for (auto &node : nodeVector) {
    auto addDofs = node->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), addDofs.begin(), addDofs.end());
  }
}


auto FaceConstraint::getTangentG1(PointerCollection &pointers, IntegrationPoint &IntegrationPt)
      -> Types::Vector3<prec>{
        return pointers.getGeometryData()->getFace(this->face)->getTangent_G1(pointers, IntegrationPt);
      }
  auto FaceConstraint::getTangentG2(PointerCollection &pointers, IntegrationPoint &IntegrationPt)
      -> Types::Vector3<prec>{
        
        return pointers.getGeometryData()->getFace(this->face)->getTangent_G2(pointers, IntegrationPt);
      }

} // namespace HierAMuS::FiniteElement
