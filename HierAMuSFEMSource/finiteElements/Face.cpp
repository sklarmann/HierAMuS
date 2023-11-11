// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "Face.h"
#include "geometry/GeometryBaseData.h"
#include "geometry/GeometryData.h"
#include "geometry/Faces/FacesRuntime.h"
#include "geometry/Faces/FacesData.h"
#include "geometry/Faces/FacesH1Interface.h"
#include "geometry/Faces/FacesHDivInterface.h"
#include "solver/GenericSolutionState.h"
#include <vector>

namespace HierAMuS::FiniteElement {
Face::~Face() = default;

auto Face::getType() -> Elementtypes { return Elementtypes::Face; }
void Face::set_pointers(PointerCollection &pointers) {
  m_face_object = pointers.getGeometryData()->getFaceData(m_face);
  m_face_runtime = pointers.getGeometryData()->getFaceRuntime(m_face);
}
void Face::setFace(indexType faceIn) { this->m_face = faceIn; }
auto Face::getVertexIds(PointerCollection &pointers) -> std::vector<indexType> {
  return m_face_runtime->getVertexNumbers();
}
auto Face::getVertex(Face::ptrCol &pointers, indexType localNumber)
    -> Geometry::VertexData & {
  std::vector<indexType> vertIds;
  vertIds = this->getVertexIds(pointers);
  
  return pointers.getGeometryData()->getVertexData(vertIds[localNumber]);
}
auto Face::getEdge(Face::ptrCol &pointers, indexType localNumber)
    -> Geometry::EdgesData & {
  auto EdgeNums = m_face_runtime->getEdgeNumbers();
  return pointers.getGeometryData()->getEdgeData(EdgeNums[localNumber]);
}
auto Face::getFace(Face::ptrCol &pointers, indexType localNumber)
    -> Geometry::FacesData * {
  return m_face_object;
}

void Face::setAllNodeBoundaryConditionMeshId(Face::ptrCol &pointers,
                                             indexType meshId, indexType dof) {
  m_face_runtime->setAllNodeBoundaryConditionMeshId(meshId, dof);
}

auto Face::getNumberOfVertices(PointerCollection &pointers) -> indexType {
  return m_face_runtime->getNumberOfVerts();
}
auto Face::getNumberOfEdges(PointerCollection &pointers) -> indexType {
  return m_face_runtime->getNumberOfEdges();
}

auto Face::getJacobian(ptrCol &pointers, IntegrationPoint &IntegrationPt)
    -> Types::Matrix22<prec> {
  return m_face_runtime->getJacobian(IntegrationPt);
}

void Face::setH1Shapes(Face::ptrCol &pointers, indexType meshid,
                       indexType order) {
  m_face_runtime->getH1Face()->setH1Shapes(meshid, order,
                                           NodeTypes::displacement);
}
void Face::getH1Dofs(Face::ptrCol &pointers,
                     std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                     indexType order) {
  Dofs = m_face_runtime->getH1Face()
             ->getH1Dofs(meshID, order);
}

auto Face::getH1Nodes(ptrCol &pointers, indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  return m_face_runtime->getH1Face()
      ->getH1Nodes(meshID, order);
}

auto Face::getH1Shapes(ptrCol &pointers, indexType order,
                       Types::Matrix22<prec> &jacobi,
                       IntegrationPoint &IntegrationPt) -> Geometry::H1Shapes {
  auto shapes = m_face_runtime->getH1Face()->getH1Shapes(order, IntegrationPt);
  Types::Matrix22<prec> jacinv = jacobi.inverse().transpose();
  shapes.shapeDeriv = jacinv * shapes.shapeDeriv;
  return shapes;
}

auto Face::getIntegrationPoints(ptrCol &pointers) -> IntegrationPoints {
  return m_face_runtime->getIntegrationPoints(this->m_id);
}

void Face::setHDivShapes(Face::ptrCol &pointers, indexType meshid,
                         indexType order, NodeTypes type) {
  m_face_runtime->getHDivFace()->setHDivShapes(meshid, order, type);
}
void Face::getHDivDofs(Face::ptrCol &pointers,
                       std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                       indexType order) {
  m_face_runtime->getHDivFace()
      ->getHDivDofs(Dofs, meshID, order, NodeTypes::displacement);
}

auto Face::getHDivShapes(PointerCollection &pointers, indexType order,
                         Types::Matrix22<prec> &jacobi,
                         IntegrationPoint &IntegrationPt)
    -> Geometry::HDivShapes {
  auto shapes = m_face_runtime->getHDivFace()->getHDivShapes(order,
                                                             IntegrationPt);
  shapes.shapes = prec(0.5) / jacobi.determinant() * jacobi * shapes.shapes;
  shapes.shapeDeriv /= jacobi.determinant() * prec(0.5);
  return shapes;
};

// HCurl Shapes

// Special Plate shapes
void Face::setSpecialPlateShapes(PointerCollection &pointers, indexType meshid,
                                 indexType order) {

  m_face_runtime->setSpecialPlateShapes(meshid, order,
                                 NodeTypes::displacement);
};
auto Face::getSpecialPlateDofs(PointerCollection &pointers, indexType meshID,
                               indexType order)
    -> std::vector<DegreeOfFreedom *> {
  std::vector<DegreeOfFreedom *> Dofs;
  Dofs = m_face_runtime->getSpecialPlateDofs(meshID, order,
                                      NodeTypes::displacement);

  return Dofs;
};
auto Face::getSpecialPlateShapes(PointerCollection &pointers, indexType order,
                                 Types::Matrix22<prec> &jacobi,
                                 IntegrationPoint &intPoint)
    -> Geometry::SpecialPlateShapes {

  // Abruf shapes in xi, eta
  auto shapes = m_face_runtime->getSpecialPlateShapes(intPoint, order);

  // TODO implemtierung der Transformation

  return shapes;
};

auto Face::getL2Shapes(ptrCol &pointers, indexType order,
                       Types::Matrix22<prec> &jacobi,
                       IntegrationPoint &IntegrationPt) -> Geometry::L2Shapes {

  return m_face_runtime->getL2Shapes(order, IntegrationPt);
}

void Face::geometryToParaview(PointerCollection &pointers,
                              vtkPlotInterface &paraviewAdapter,
                              indexType mainMesh, indexType subMesh) {
  m_face_runtime->geometryToParaview(paraviewAdapter, mainMesh, subMesh);
};

void Face::computeWeightsParaview(PointerCollection &pointers,
                                  vtkPlotInterface &paraviewAdapter,
                                  indexType mainMesh, indexType subMesh) {
  m_face_runtime->computeWeightsParaview(paraviewAdapter, mainMesh, subMesh);
}

void Face::H1SolutionToParaview(PointerCollection &pointers,
                                vtkPlotInterface &paraviewAdapter,
                                indexType mainMesh, indexType subMesh,
                                indexType meshId, indexType order,
                                std::string name) {
  auto Dofs = m_face_runtime->getH1Face()->getH1Dofs(meshId, order);
  Types::VectorX<prec> sol = pointers.getSolutionState()->getSolution(Dofs);
  m_face_runtime->H1SolutionToParaview(paraviewAdapter, mainMesh, subMesh,
                                order, sol, name);
}

void Face::projectDataToParaviewVertices(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {
  m_face_runtime->projectDataToParaviewVertices(paraviewAdapter, mainMesh,
                                         subMesh, order, IntegrationPt, data,
                                         numberComponents, name);
}

void Face::getL2Dofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                     indexType meshID, indexType order) {
  m_face_runtime->getL2Dofs(Dofs, meshID, order);
}

void Face::setL2Shapes(ptrCol &pointers, indexType meshid, indexType order) {
  m_face_runtime->setL2Shapes(meshid, order, NodeTypes::displacement);
}

} // namespace HierAMuS::FiniteElement
