// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "plot/vtkplotClassBase.h"


#include <finiteElements/Edge.h>
#include <finiteElements/GenericFiniteElement.h>

#include <pointercollection/pointercollection.h>

#include <geometry/Base.h>
#include <geometry/Edges.h>
#include <geometry/Vertex.h>

#include <geometry/GeometryData.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>

#include <equations/NodeSet.h>
#include <shapefunctions/LagrangeShape.h>
#include <solver/GenericSolutionState.h>
#include <solver/SolutionTypes.h>

namespace HierAMuS::FiniteElement {
  

Edge::~Edge() = default;

void Edge::setAllNodeBoundaryConditionMeshId(ptrCol &pointers, indexType meshId,
                                             indexType dof) {

  auto &tempEdge = pointers.getGeometryData()->getEdge(m_edge);
  tempEdge.setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);

  indexType numVerts = tempEdge.getNumberOfVerts();
  for (auto i = 0; i < numVerts; ++i) {
    auto &Vert = tempEdge.getVertex(pointers, i);
    Vert.setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);
  }
}

auto Edge::getA1Vector(ptrCol &pointers, IntegrationPoint &integration_point)
    -> Types::Vector3<prec> {
  auto &edge = pointers.getGeometryData()->getEdge(m_edge);
  return edge.getA1Vector(pointers,integration_point);
}

void Edge::setH1Shapes(ptrCol &pointers, indexType meshID, indexType order) {

  auto &tempEdge = pointers.getGeometryData()->getEdge(m_edge);
  tempEdge.setH1Shapes(pointers, meshID, order, NodeTypes::displacement);
}

void Edge::setEdges(std::vector<indexType> &edgesIn) {
  if (edgesIn.size() == 1) {
    m_edge = edgesIn[0];
  }
}

auto Edge::getVertexId(ptrCol &pointers, indexType num) -> indexType {
  auto geom = pointers.getGeometryData();
  auto &geomEle = geom->getEdge(m_edge);
  std::vector<indexType> verts;
  geomEle.getVerts(verts);

  return verts[num];
}

auto Edge::getVertex(ptrCol &pointers, indexType num) -> Geometry::Vertex & {
  auto geodata = pointers.getGeometryData();
  auto &geotemp = geodata->getEdge(m_edge);
  std::vector<indexType> verts;
  geotemp.getVerts(verts);
  if (num >= verts.size()) {
    throw std::runtime_error("Error in Edge requested Vertex not in list!");
  }
  return geodata->getVertex(verts[num]);
}

auto Edge::getEdge(Edge::ptrCol &pointers, indexType localNumber)
    -> Geometry::Edges & {
  return pointers.getGeometryData()->getEdge(m_edge);
}


void Edge::getH1Dofs(ptrCol &pointers,
                            std::vector<DegreeOfFreedom *> &Dofs,
                            indexType meshID, indexType order) {
  auto &edge = pointers.getGeometryData()->getEdge(m_edge);
  edge.getH1Dofs(pointers, Dofs, meshID, order);
}

auto Edge::getIntegrationPoints(ptrCol &pointers)
    -> IntegrationPoints {
  auto &edge = pointers.getGeometryData()->getEdge(m_edge);
  auto intPoints = edge.getIntegrationPoints(pointers,this->id);
  return intPoints;
}

void Edge::getJacobian(ptrCol &pointers, prec &jacobi, prec xsi) {
  auto &edge = pointers.getGeometryData()->getEdge(m_edge);
  jacobi = edge.getJacobian(pointers, xsi);
}

auto Edge::getJacobian(ptrCol &pointers, IntegrationPoint &IntegrationPt)
    -> Types::MatrixXX<prec> {
  auto &edge = pointers.getGeometryData()->getEdge(m_edge);

  return edge.getJacobian(pointers,IntegrationPt);
}

void Edge::getH1Shapes(ptrCol &pointers, indexType order, prec jacobi,
                              Types::VectorX<prec> &shape,
                              Types::VectorX<prec> &shapeDerivative, prec xsi) {
  auto &edge = pointers.getGeometryData()->getEdge(m_edge);
  IntegrationPoint integration_point;
  integration_point.xi = xsi;
  auto shapes = edge.getH1Shapes(pointers, order, integration_point);
  shape = shapes.shapes;
  shapeDerivative = shapes.shapeDeriv / jacobi;
}

auto Edge::getH1Shapes(ptrCol &pointers, indexType order,
                       Types::MatrixXX<prec> &jacobi,
                       IntegrationPoint &IntegrationPt) -> Geometry::H1Shapes {

  auto &edge = pointers.getGeometryData()->getEdge(m_edge);
  auto shapes = edge.getH1Shapes(pointers, order, IntegrationPt);
  shapes.shapeDeriv /= jacobi(0, 0);

  return shapes;
}


void Edge::geometryToParaview(PointerCollection &pointers,
                              vtkPlotInterface &paraviewAdapter,
                              indexType mainMesh, indexType subMesh) {

  auto &temp = pointers.getGeometryData()->getEdge(this->m_edge);
  temp.geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
}
void Edge::computeWeightsParaview(PointerCollection &pointers,
                                  vtkPlotInterface &paraviewAdapter,
                                  indexType mainMesh, indexType subMesh) {
  auto &temp = pointers.getGeometryData()->getEdge(this->m_edge);
  temp.geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
}
void Edge::H1SolutionToParaview(PointerCollection &pointers,
                                vtkPlotInterface &paraviewAdapter,
                                indexType mainMesh, indexType subMesh,
                                indexType meshId, indexType order,
                                std::string name)
{
  auto &geoElem = pointers.getGeometryData()->getEdge(this->m_edge);
  geoElem.H1SolutionToParaview(pointers, paraviewAdapter, mainMesh, subMesh,
                                meshId, order, name);
}
void Edge::projectDataToParaviewVertices(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name)
{

  auto &geoElem = pointers.getGeometryData()->getEdge(this->m_edge);
  geoElem.projectDataToParaviewVertices(pointers, paraviewAdapter, mainMesh,
                                         subMesh, order, IntegrationPt, data,
                                         numberComponents, name);
};

} // namespace HierAMuS
