// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "plot/vtkplotClassBase.h"

#include <finiteElements/Edge.h>
#include <finiteElements/GenericFiniteElement.h>

#include <pointercollection/pointercollection.h>

#include "geometry/Edges/EdgeH1ShapesInterface.h"
#include "geometry/Edges/EdgesRuntime.h"
#include "geometry/VertexRuntime.h"
#include <geometry/Edges/EdgesData.h>
#include <geometry/GeometryBaseData.h>

#include <geometry/GeometryData.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>

#include <shapefunctions/LagrangeShape.h>
#include <solver/GenericSolutionState.h>
#include <solver/SolutionTypes.h>

namespace HierAMuS::FiniteElement {

Edge::~Edge() = default;

void Edge::setAllNodeBoundaryConditionMeshId(ptrCol &pointers, indexType meshId,
                                             indexType dof) {

  m_edge_runtime->setAllNodeBoundaryConditionMeshId(meshId, dof);

  indexType numVerts = m_edge_runtime->getNumberOfVerts();
  for (auto i = 0; i < numVerts; ++i) {
    auto &Vert = m_edge_runtime->getVertex(i);
    Vert.setAllNodeBoundaryConditionMeshId(meshId, dof);
  }
}

auto Edge::getA1Vector(ptrCol &pointers, IntegrationPoint &integration_point)
    -> Types::Vector3<prec> {
  return m_edge_runtime->getA1Vector(integration_point);
}

void Edge::setH1Shapes(ptrCol &pointers, indexType meshID, indexType order) {

  m_edge_runtime->getH1Edge()->setH1Shapes(meshID, order,
                                           NodeTypes::displacement);
}

void Edge::setEdges(std::vector<indexType> &edgesIn) {
  if (edgesIn.size() == 1) {
    m_edge = edgesIn[0];
  }
}

auto Edge::getVertexId(ptrCol &pointers, indexType num) -> indexType {
  auto verts = m_edge_runtime->getVertexNumbers();
  return verts[num];
}

auto Edge::getVertex(ptrCol &pointers, indexType num)
    -> Geometry::VertexData & {
  auto geodata = pointers.getGeometryData();
  auto verts = m_edge_runtime->getVertexNumbers();
  if (num >= static_cast<indexType>(verts.size())) {
    throw std::runtime_error("Error in Edge requested Vertex not in list!");
  }
  return geodata->getVertexData(verts[num]);
}

auto Edge::getEdge(Edge::ptrCol &pointers, indexType localNumber)
    -> Geometry::EdgesData & {
  return pointers.getGeometryData()->getEdgeData(m_edge);
}

void Edge::getH1Dofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                     indexType meshID, indexType order) {
  m_edge_runtime->getH1Edge()->getH1Dofs(Dofs, meshID, order);
}

auto Edge::getIntegrationPoints(ptrCol &pointers) -> IntegrationPoints {
  auto &edge = pointers.getGeometryData()->getEdgeData(m_edge);
  auto intPoints = edge.getIntegrationPoints(this->m_id);
  return intPoints;
}

auto Edge::getJacobian(ptrCol &pointers, IntegrationPoint &IntegrationPt)
    -> prec {
  auto &edge = pointers.getGeometryData()->getEdgeData(m_edge);
  auto jj = edge.getJacobian(IntegrationPt);
  return jj;
}

auto Edge::getH1Shapes(ptrCol &pointers, indexType order, prec jacobi,
                       IntegrationPoint &IntegrationPt) -> Geometry::H1Shapes {

  auto shapes = m_edge_runtime->getH1Edge()->getH1Shapes(order, IntegrationPt);
  shapes.shapeDeriv /= jacobi;

  return shapes;
}

void Edge::geometryToParaview(PointerCollection &pointers,
                              vtkPlotInterface &paraviewAdapter,
                              indexType mainMesh, indexType subMesh) {

  auto &temp = *pointers.getGeometryData()->getEdgeRuntime(this->m_edge);
  temp.geometryToParaview(paraviewAdapter, mainMesh, subMesh);
}
void Edge::computeWeightsParaview(PointerCollection &pointers,
                                  vtkPlotInterface &paraviewAdapter,
                                  indexType mainMesh, indexType subMesh) {
  auto &temp = *pointers.getGeometryData()->getEdgeRuntime(this->m_edge);
  temp.geometryToParaview(paraviewAdapter, mainMesh, subMesh);
}
void Edge::H1SolutionToParaview(PointerCollection &pointers,
                                vtkPlotInterface &paraviewAdapter,
                                indexType mainMesh, indexType subMesh,
                                indexType meshId, indexType order,
                                std::string name) {

  auto &geoElem = *pointers.getGeometryData()->getEdgeRuntime(this->m_edge);
  std::vector<DegreeOfFreedom *> Dofs;
  geoElem.getH1Edge()->getH1Dofs(Dofs, meshId, order);
  Types::VectorX<prec> sol = pointers.getSolutionState()->getSolution(Dofs);

  geoElem.H1SolutionToParaview(paraviewAdapter, mainMesh, subMesh, order, sol,
                               name);
}
void Edge::projectDataToParaviewVertices(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {

  auto &geoElem = *pointers.getGeometryData()->getEdgeRuntime(this->m_edge);
  geoElem.projectDataToParaviewVertices(paraviewAdapter, mainMesh, subMesh,
                                        order, IntegrationPt, data,
                                        numberComponents, name);
};

void Edge::set_pointers(PointerCollection &pointers) {
  m_edge_pointer = &pointers.getGeometryData()->getEdgeData(m_edge);
  m_edge_runtime = pointers.getGeometryData()->getEdgeRuntime(m_edge);
}

auto Edge::getElementType() -> Elementtypes { return Elementtypes::Edge; }

} // namespace HierAMuS::FiniteElement
