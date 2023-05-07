// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "Face.h"
#include "equations/Nodetypes.h"
#include "geometry/Base.h"
#include <vector>

namespace HierAMuS::FiniteElement {
Face::~Face() = default;

auto Face::getType() -> Elementtypes { return Elementtypes::Face; }
void Face::setFace(indexType faceIn) { this->face = faceIn; }
auto Face::getVertexIds(PointerCollection& pointers) -> std::vector<indexType> {
  auto tempFace =
      pointers.getGeometryData()->getFace(this->face);
  std::vector<indexType> vertIds;
  tempFace->getVerts(vertIds);
  return vertIds;
}
auto Face::getVertex(Face::ptrCol &pointers, indexType localNumber)
    -> Geometry::Vertex & {
  std::vector<indexType> vertIds;
  vertIds = this->getVertexIds(pointers);

  return pointers.getGeometryData()->getVertex(vertIds[localNumber]);
}
auto Face::getEdge(Face::ptrCol &pointers, indexType localNumber)
    -> Geometry::Edges & {
  Geometry::Base *temp;
  temp = pointers.getGeometryData()->getFace(this->face);
  std::vector<indexType> EdgeNums;
  temp->getEdges(EdgeNums);
  return pointers.getGeometryData()->getEdge(EdgeNums[localNumber]);
}
auto Face::getFace(Face::ptrCol &pointers, indexType localNumber)
    -> Geometry::Faces * {
  return pointers.getGeometryData()->getFace(this->face);
}

void Face::setAllNodeBoundaryConditionMeshId(Face::ptrCol &pointers,
                                             indexType meshId, indexType dof) {
  auto temp = pointers.getGeometryData()->getFace(this->face);
  temp->setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);
  // GenericFiniteElement::setAllNodeBoundaryConditionMeshId(pointers, meshId,
  // dof);
}

auto Face::getNumberOfVertices(PointerCollection& pointers) -> indexType {
  return GenericFiniteElement::getNumberOfVertices(pointers);
}
auto Face::getNumberOfEdges(PointerCollection& pointers) -> indexType {
  return GenericFiniteElement::getNumberOfEdges(pointers);
}

auto Face::getJacobian(ptrCol &pointers, IntegrationPoint &IntegrationPt)
    -> Types::MatrixXX<prec> {
  auto geomElement = pointers.getGeometryData()->getFace(this->face);
  return geomElement->getJacobian(pointers, IntegrationPt);
}
void Face::getJacobian(Face::ptrCol &pointers, Types::Matrix22<prec> &jacobi,
                       prec xsi, prec eta) {
  GenericFiniteElement::getJacobian(pointers, jacobi, xsi, eta);
}
void Face::setH1Shapes(Face::ptrCol &pointers, indexType meshid,
                       indexType order) {
  Geometry::Faces *tempFace;
  tempFace = pointers.getGeometryData()->getFace(this->face);
  tempFace->setH1Shapes(pointers, meshid, order, NodeTypes::displacement);
}
void Face::getH1Dofs(Face::ptrCol &pointers,
                     std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                     indexType order) {
  pointers.getGeometryData()
      ->getFace(this->face)
      ->getH1Dofs(pointers, Dofs, meshID, order);
}

auto Face::getH1Nodes(ptrCol &pointers, indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  return pointers.getGeometryData()
      ->getFace(this->face)
      ->getH1Nodes(pointers, meshID, order);
}

void Face::getH1Shapes(ptrCol &pointers, indexType order,
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
void Face::getH1Shapes(Face::ptrCol &pointers, indexType order,
                       Types::Matrix22<prec> &jacobi,
                       Types::VectorX<prec> &shape,
                       Types::Matrix2X<prec> &dshape, prec xi, prec eta) {
  auto geoFace = pointers.getGeometryData()->getFace(this->face);
  geoFace->getH1Shapes(pointers, order, shape, dshape, xi, eta);
  dshape = jacobi.inverse().transpose() * dshape;
}

auto Face::getH1Shapes(ptrCol &pointers, indexType order,
                       Types::MatrixXX<prec> &jacobi,
                       IntegrationPoint &IntegrationPt) -> Geometry::H1Shapes {
  auto geoFace = pointers.getGeometryData()->getFace(this->face);
  auto shapes = geoFace->getH1Shapes(pointers, order, IntegrationPt);
  shapes.shapeDeriv = jacobi.inverse().transpose() * shapes.shapeDeriv;
  return shapes;
}

auto Face::getIntegrationPoints(ptrCol &pointers) -> IntegrationPoints {
  auto temp = pointers.getGeometryData()->getFace(this->face);
  return temp->getIntegrationPoints(pointers,this->id);
}

void Face::setHDivShapes(Face::ptrCol &pointers, indexType meshid,
                         indexType order, NodeTypes type) {
  Geometry::Faces *tempFace;
  tempFace = pointers.getGeometryData()->getFace(this->face);
  tempFace->setHDivShapes(pointers, meshid, order, type);
}
void Face::getHDivDofs(Face::ptrCol &pointers,
                       std::vector<DegreeOfFreedom *> &Dofs,
                       indexType meshID, indexType order) {
  pointers.getGeometryData()
      ->getFace(this->face)
      ->getHDivDofs(pointers, Dofs, meshID, order, NodeTypes::displacement);
}
void Face::getHDivShapes(Face::ptrCol &pointers, indexType order,
                         Types::Matrix22<prec> &jacobi,
                         Types::Matrix2X<prec> &shape,
                         Types::VectorX<prec> &dshape, prec xi,
                         prec eta) {
  Geometry::Faces *tempFace;
  tempFace = pointers.getGeometryData()->getFace(this->face);
  tempFace->getHDivShapes(pointers, order, shape, dshape, xi, eta);
  shape = prec(0.5) / jacobi.determinant() * jacobi * shape;
  dshape = dshape / jacobi.determinant() * prec(0.5);
}

auto Face::getHDivShapes(PointerCollection &pointers, indexType order,
                         Types::MatrixXX<prec> &jacobi,
                         IntegrationPoint &IntegrationPt)
    -> Geometry::HDivShapes {
  auto geoFace = pointers.getGeometryData()->getFace(this->face);
  auto shapes = geoFace->getHDivShapes(pointers, order, IntegrationPt);
  shapes.shapes =
      prec(0.5) / jacobi.determinant() * jacobi * shapes.shapes;
  shapes.shapeDeriv /=  jacobi.determinant() * prec(0.5);
  return shapes;
};

// HCurl Shapes

// Special Plate shapes
void Face::setSpecialPlateShapes(PointerCollection &pointers, indexType meshid,
                                 indexType order) {

  auto geoFace = pointers.getGeometryData()->getFace(this->face);
  geoFace->setSpecialPlateShapes(pointers, meshid, order,
                                 NodeTypes::displacement);
};
auto Face::getSpecialPlateDofs(PointerCollection &pointers, indexType meshID,
                               indexType order)
    -> std::vector<DegreeOfFreedom *> {
  std::vector<DegreeOfFreedom *> Dofs;
  auto geoFace = pointers.getGeometryData()->getFace(this->face);
  Dofs = geoFace->getSpecialPlateDofs(pointers, meshID, order,
                                      NodeTypes::displacement);

  return Dofs;
};
auto Face::getSpecialPlateShapes(PointerCollection &pointers, indexType order,
                                 Types::MatrixXX<prec> &jacobi,
                                 IntegrationPoint &intPoint)
    -> Geometry::SpecialPlateShapes {
  // Geometry::SpecialPlateShapes shapes;
  auto geoFace = pointers.getGeometryData()->getFace(this->face);

  // Abruf shapes in xi, eta
  auto shapes = geoFace->getSpecialPlateShapes(pointers, intPoint, order);

  // TODO implemtierung der Transformation

  return shapes;
};
void Face::toParaviewAdapter(PointerCollection &ptrCol,
                             vtkPlotInterface &catalyst,
                             const ParaviewSwitch &ToDo) {}

void Face::geometryToParaview(PointerCollection &pointers,
                              vtkPlotInterface &paraviewAdapter,
                              indexType mainMesh, indexType subMesh) {
  auto temp = pointers.getGeometryData()->getFace(this->face);
  temp->geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
};

void Face::computeWeightsParaview(PointerCollection &pointers,
                                  vtkPlotInterface &paraviewAdapter,
                                  indexType mainMesh, indexType subMesh) {
  auto temp = pointers.getGeometryData()->getFace(this->face);
  temp->computeWeightsParaview(pointers, paraviewAdapter, mainMesh, subMesh);
}

void Face::H1SolutionToParaview(PointerCollection &pointers,
                                vtkPlotInterface &paraviewAdapter,
                                indexType mainMesh, indexType subMesh,
                                indexType meshId, indexType order,
                                std::string name) {
  auto geoElem = pointers.getGeometryData()->getFace(this->face);
  geoElem->H1SolutionToParaview(pointers, paraviewAdapter, mainMesh, subMesh,
                                meshId, order, name);
}

void Face::projectDataToParaviewVertices(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {
  auto geoElem = pointers.getGeometryData()->getFace(this->face);
  geoElem->projectDataToParaviewVertices(pointers, paraviewAdapter, mainMesh,
                                         subMesh, order, IntegrationPt, data,
                                         numberComponents, name);
}

} // namespace HierAMuS::FiniteElement
