// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "Volume.h"
#include "MatrixTypes.h"
#include "geometry/Base.h"
#include <vector>

namespace HierAMuS::FiniteElement {
Volume::~Volume() = default;

auto Volume::getType() -> Elementtypes { return Elementtypes::Volume; }
void Volume::setVolume(indexType volumeIn) { this->volume = volumeIn; }
auto Volume::getVertexIds(PointerCollection& pointers) -> std::vector<indexType> {
  auto tempVol =
      pointers.getGeometryData()->getVolume(this->volume);
  std::vector<indexType> vertIds;
  tempVol->getVerts(vertIds);
  return vertIds;
}

auto Volume::getVertex(Volume::ptrCol &pointers, indexType localNumber)
    -> Geometry::Vertex & {
  std::vector<indexType> vertIds;
  vertIds = this->getVertexIds(pointers);

  return pointers.getGeometryData()->getVertex(vertIds[localNumber]);
}

auto Volume::getEdge(Volume::ptrCol &pointers, indexType localNumber)
    -> Geometry::Edges & {
  auto temp = pointers.getGeometryData()->getVolume(this->volume);
  std::vector<indexType> EdgeNums;
  temp->getEdges(EdgeNums);
  return pointers.getGeometryData()->getEdge(EdgeNums[localNumber]);
}

auto Volume::getFace(Volume::ptrCol &pointers, indexType localNumber)
    -> Geometry::Faces * {
  auto temp = pointers.getGeometryData()->getVolume(this->volume);
  std::vector<indexType> FaceNums;
  temp->getFaces(FaceNums);
  return pointers.getGeometryData()->getFace(FaceNums[localNumber]);
}


auto Volume::getVolume(ptrCol &pointers, indexType localNumber) -> Geometry::Volumes *
{
  return pointers.getGeometryData()->getVolume(this->volume);
}



void Volume::setAllNodeBoundaryConditionMeshId(Volume::ptrCol &pointers,
                                             indexType meshId, indexType dof) {
  auto temp = pointers.getGeometryData()->getVolume(this->volume);
  temp->setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);
}

auto Volume::getNumberOfVertices(PointerCollection& pointers) -> indexType {
  auto temp = pointers.getGeometryData()->getVolume(this->volume);
  return temp->getNumberOfVerts();
}

auto Volume::getNumberOfEdges(PointerCollection& pointers) -> indexType {
  auto temp = pointers.getGeometryData()->getVolume(this->volume);
  return temp->getNumberOfEdges();
}

auto Volume::getNumberOfFaces(PointerCollection& pointers) -> indexType {
  auto temp = pointers.getGeometryData()->getVolume(this->volume);
  return temp->getNumberOfFaces();
}

auto Volume::getJacobian(ptrCol &pointers, IntegrationPoint &IntegrationPt)
    -> Types::MatrixXX<prec> {
  auto geomElement = pointers.getGeometryData()->getVolume(this->volume);
  return geomElement->getJacobian(pointers, IntegrationPt);
}

void Volume::getJacobian(Volume::ptrCol &pointers, Types::Matrix33<prec> &jacobi,
                       prec xsi, prec eta, prec zeta) {
  auto geomElement = pointers.getGeometryData()->getVolume(this->volume);
  geomElement->getJacobian(pointers, jacobi, xsi, eta, zeta);

}

void Volume::setH1Shapes(Volume::ptrCol &pointers, indexType meshid,
                       indexType order) {

  auto tempFace = pointers.getGeometryData()->getVolume(this->volume);
  tempFace->setH1Shapes(pointers, meshid, order, NodeTypes::displacement);
}

void Volume::getH1Dofs(Volume::ptrCol &pointers,
                     std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                     indexType order) {
  pointers.getGeometryData()
      ->getVolume(this->volume)
      ->getH1Dofs(pointers, Dofs, meshID, order);
}

auto Volume::getH1Nodes(ptrCol &pointers, indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  return pointers.getGeometryData()
      ->getVolume(this->volume)
      ->getH1Nodes(pointers, meshID, order);
}

void Volume::getH1Shapes(ptrCol &pointers, indexType order,
                       Types::MatrixXX<prec> &jacobi,
                       Types::VectorX<prec> &shape,
                       Types::MatrixXX<prec> &shapeDerivative,
                       IntegrationPoint &IntegrationPt) {
  auto geomElement = pointers.getGeometryData()->getVolume(this->volume);
  auto shapes = geomElement->getH1Shapes(pointers, order, IntegrationPt);
  shape = shapes.shapes;
  shapeDerivative = shapes.shapeDeriv;
  shapeDerivative = jacobi.inverse().transpose() * shapeDerivative;
}

void Volume::getH1Shapes(Volume::ptrCol &pointers, indexType order,
                       Types::Matrix33<prec> &jacobi,
                       Types::VectorX<prec> &shape,
                       Types::Matrix3X<prec> &dshape, prec xi, prec eta,prec zeta) {
  auto geomElement = pointers.getGeometryData()->getVolume(this->volume);
  geomElement->getH1Shapes(pointers, order, shape, dshape, xi, eta, zeta);
  dshape = jacobi.inverse().transpose() * dshape;
}

auto Volume::getH1Shapes(ptrCol &pointers, indexType order,
                       Types::MatrixXX<prec> &jacobi,
                       IntegrationPoint &IntegrationPt) -> Geometry::H1Shapes {
  auto geoVol = pointers.getGeometryData()->getVolume(this->volume);
  auto shapes = geoVol->getH1Shapes(pointers, order, IntegrationPt);
  Types::Matrix33<prec> jacobiInv = jacobi;
  jacobiInv = jacobiInv.inverse().transpose();

  shapes.shapeDeriv = jacobiInv * shapes.shapeDeriv;
  return shapes;
}

auto Volume::getIntegrationPoints(ptrCol &pointers) -> IntegrationPoints {
  auto temp = pointers.getGeometryData()->getVolume(this->volume);
  return temp->getIntegrationPoints(pointers,this->id);
}

void Volume::setHDivShapes(Volume::ptrCol &pointers, indexType meshid,
                         indexType order, NodeTypes type) {
  GenericFiniteElement::setHDivShapes(pointers, meshid, order, type);
}
void Volume::getHDivDofs(Volume::ptrCol &pointers,
                       std::vector<DegreeOfFreedom *> &Dofs,
                       indexType meshID, indexType order) {
  GenericFiniteElement::getHDivDofs(pointers, Dofs, meshID, order);
}
void Volume::getHDivShapes(Volume::ptrCol &pointers, indexType order,
                           Types::Matrix22<prec> &jacobi,
                           Types::Matrix2X<prec> &shape,
                           Types::VectorX<prec> &dshape, prec xi,
                           prec eta) {
  GenericFiniteElement::getHDivShapes(pointers, order, jacobi, shape, dshape,
                                      xi, eta);
}
void Volume::toParaviewAdapter(PointerCollection &ptrCol,
                             vtkPlotInterface &catalyst,
                             const ParaviewSwitch &ToDo) {}

void Volume::geometryToParaview(PointerCollection &pointers,
                              vtkPlotInterface &paraviewAdapter,
                              indexType mainMesh, indexType subMesh) {
  auto temp = pointers.getGeometryData()->getVolume(this->volume);
  temp->geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
};

void Volume::computeWeightsParaview(PointerCollection &pointers,
                                  vtkPlotInterface &paraviewAdapter,
                                  indexType mainMesh, indexType subMesh) {
  auto temp = pointers.getGeometryData()->getVolume(this->volume);
  temp->computeWeightsParaview(pointers, paraviewAdapter, mainMesh, subMesh);
}

void Volume::H1SolutionToParaview(PointerCollection &pointers,
                                vtkPlotInterface &paraviewAdapter,
                                indexType mainMesh, indexType subMesh,
                                indexType meshId, indexType order,
                                std::string name) {
  auto geoElem = pointers.getGeometryData()->getVolume(this->volume);
  geoElem->H1SolutionToParaview(pointers, paraviewAdapter, mainMesh, subMesh,
                                meshId, order, name);
}

void Volume::projectDataToParaviewVertices(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {
  auto geoElem = pointers.getGeometryData()->getVolume(this->volume);
  geoElem->projectDataToParaviewVertices(pointers, paraviewAdapter, mainMesh,
                                         subMesh, order, IntegrationPt, data,
                                         numberComponents, name);
}

} // namespace HierAMuS
