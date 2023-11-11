// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "Volume.h"
#include "MatrixTypes.h"
#include "geometry/GeometryBaseData.h"
#include "geometry/GeometryData.h"
#include "geometry/Volumes/VolumesH1Interface.h"
#include "geometry/Volumes/VolumesRuntime.h"
#include "solver/GenericSolutionState.h"
#include <vector>

namespace HierAMuS::FiniteElement {
Volume::~Volume() = default;

auto Volume::getType() -> Elementtypes { return Elementtypes::Volume; }
void Volume::set_pointers(PointerCollection &pointers) {
  m_volume_object = pointers.getGeometryData()->getVolumeData(m_volume);
  m_volume_runtime =
      pointers.getGeometryData()->getVolumeRuntime(m_volume);
}
void Volume::setVolume(indexType volumeIn) { this->m_volume = volumeIn; }
auto Volume::getVertexIds(PointerCollection &pointers)
    -> std::vector<indexType> {
  return m_volume_object->getVertexNumbers();
}

auto Volume::getVertex(Volume::ptrCol &pointers, indexType localNumber)
    -> Geometry::VertexData & {
  std::vector<indexType> vertIds;
  vertIds = this->getVertexIds(pointers);

  return pointers.getGeometryData()->getVertexData(vertIds[localNumber]);
}

auto Volume::getEdge(Volume::ptrCol &pointers, indexType localNumber)
    -> Geometry::EdgesData & {
  return *m_volume_object->getEdge(localNumber);
}

auto Volume::getFace(Volume::ptrCol &pointers, indexType localNumber)
    -> Geometry::FacesData * {
  std::vector<indexType> FaceNums;
  m_volume_object->getFaceNumbers(FaceNums);
  return pointers.getGeometryData()->getFaceData(FaceNums[localNumber]);
}

auto Volume::getVolume(ptrCol &pointers, indexType localNumber)
    -> Geometry::VolumesData * {
  return m_volume_object;
}

void Volume::setAllNodeBoundaryConditionMeshId(Volume::ptrCol &pointers,
                                               indexType meshId,
                                               indexType dof) {
  m_volume_object->setAllNodeBoundaryConditionMeshId(meshId, dof);
}

auto Volume::getNumberOfVertices(PointerCollection &pointers) -> indexType {
  return m_volume_runtime->getNumberOfVerts();
}

auto Volume::getNumberOfEdges(PointerCollection &pointers) -> indexType {
  return m_volume_runtime->getNumberOfEdges();
}

auto Volume::getNumberOfFaces(PointerCollection &pointers) -> indexType {
  return m_volume_runtime->getNumberOfFaces();
}

auto Volume::getJacobian(ptrCol &pointers, IntegrationPoint &IntegrationPt)
    -> Types::Matrix33<prec> {
  return m_volume_runtime->getJacobian(IntegrationPt);
}

void Volume::setH1Shapes(Volume::ptrCol &pointers, indexType meshid,
                         indexType order) {

  m_volume_object->setH1Shapes(meshid, order, NodeTypes::displacement);
}

void Volume::getH1Dofs(Volume::ptrCol &pointers,
                       std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                       indexType order) {
  m_volume_runtime->getH1Volume()->getH1Dofs(Dofs, meshID, order);
}

auto Volume::getH1Nodes(ptrCol &pointers, indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  return m_volume_runtime->getH1Volume()->getH1Nodes(meshID, order);
}

auto Volume::getH1Shapes(ptrCol &pointers, indexType order,
                         Types::Matrix33<prec> &jacobi,
                         IntegrationPoint &IntegrationPt)
    -> Geometry::H1Shapes {
  auto shapes =
      m_volume_runtime->getH1Volume()->getH1Shapes(order, IntegrationPt);
  Types::Matrix33<prec> jacobiInv = jacobi.inverse().transpose();

  shapes.shapeDeriv = jacobiInv * shapes.shapeDeriv;
  return shapes;
}

auto Volume::getIntegrationPoints(ptrCol &pointers) -> IntegrationPoints {
  return m_volume_runtime->getIntegrationPoints(this->m_id);
}

void Volume::setHDivShapes(Volume::ptrCol &pointers, indexType meshid,
                           indexType order, NodeTypes type) {
  GenericFiniteElement::setHDivShapes(pointers, meshid, order, type);
}
void Volume::getHDivDofs(Volume::ptrCol &pointers,
                         std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) {
  GenericFiniteElement::getHDivDofs(pointers, Dofs, meshID, order);
}

void Volume::geometryToParaview(PointerCollection &pointers,
                                vtkPlotInterface &paraviewAdapter,
                                indexType mainMesh, indexType subMesh) {
  m_volume_runtime->geometryToParaview(paraviewAdapter, mainMesh,
                                       subMesh);
};

void Volume::computeWeightsParaview(PointerCollection &pointers,
                                    vtkPlotInterface &paraviewAdapter,
                                    indexType mainMesh, indexType subMesh) {
  m_volume_runtime->computeWeightsParaview(paraviewAdapter, mainMesh,
                                           subMesh);
}

void Volume::H1SolutionToParaview(PointerCollection &pointers,
                                  vtkPlotInterface &paraviewAdapter,
                                  indexType mainMesh, indexType subMesh,
                                  indexType meshId, indexType order,
                                  std::string name) {
  std::vector<DegreeOfFreedom *> Dofs;
  m_volume_runtime->getH1Volume()->getH1Dofs(Dofs, meshId, order);
  Types::VectorX<prec> sol = pointers.getSolutionState()->getSolution(Dofs);
  m_volume_runtime->H1SolutionToParaview(paraviewAdapter, mainMesh,
                                         subMesh, order, sol, name);
}

void Volume::projectDataToParaviewVertices(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {
  m_volume_runtime->projectDataToParaviewVertices(
      paraviewAdapter, mainMesh, subMesh, order, IntegrationPt, data,
      numberComponents, name);
}

} // namespace HierAMuS::FiniteElement
