// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "GenericNodes.h"

#include "geometry/VertexData.h"
#include <geometry/GeometryTypes.h>
#include <geometry/VertexRuntime.h>


#include "LoadList.h"

#include <vector>

#include <plot/vtkplotClass.h>

namespace HierAMuS::Geometry {

VertexRuntime::VertexRuntime(VertexData &data)
    : m_Vertex_Data_Element(&data), GeometryBaseRuntime(data),
      coor(data.getCoordinates()) {}

VertexRuntime::VertexRuntime(VertexRuntime &other)
    : GeometryBaseRuntime(other),
      m_Vertex_Data_Element(other.m_Vertex_Data_Element), coor(other.coor) {}

VertexRuntime::VertexRuntime(VertexRuntime &&other)
    : GeometryBaseRuntime(other),
      m_Vertex_Data_Element(other.m_Vertex_Data_Element),
      coor(std::move(other.coor)) {}

VertexRuntime &VertexRuntime::operator=(VertexRuntime &other) {
  GeometryBaseRuntime::operator=(other);
  coor = other.coor;
  m_Vertex_Data_Element = other.m_Vertex_Data_Element;
  return *this;
}

VertexRuntime &VertexRuntime::operator=(VertexRuntime &&other) {
  GeometryBaseRuntime::operator=(other);
  m_Vertex_Data_Element = other.m_Vertex_Data_Element;
  coor = std::move(other.coor);
  return *this;

}



VertexRuntime::~VertexRuntime() {}

auto VertexRuntime::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::VertexRuntime::type;
}

void VertexRuntime::setCoordinates(prec x, prec y, prec z) {
  m_Vertex_Data_Element->setCoordinates(x, y, z);
}

void VertexRuntime::setCoordinates(const Types::Vector3<prec> &coorIn) {
  m_Vertex_Data_Element->setCoordinates(coorIn);
}

void VertexRuntime::print(spdlog::logger &Log) {
  Log.debug(*this);
  this->printEqInfo(Log);
}

auto VertexRuntime::getCoordinates() -> Types::Vector3<prec> {
  return m_Vertex_Data_Element->getCoordinates();
}

void VertexRuntime::setBoundaryCondition(indexType meshId, indexType order,
                                         ShapeFunctionTypes shapeType,
                                         Types::Vector3<indexType> &dofs,
                                         bool set) {

  std::vector<GenericNodes *> tnodes;
  this->getNodes(tnodes, meshId);

  if (set) {
    GenericNodes *tnode;
    tnode = tnodes[0];
    for (auto i = 0; i < 3; ++i) {
      if (dofs(i) != 0) {
        tnode->setBoundaryCondition(i);
      } else {
        tnode->unsetBoundaryCondition(i);
      }
    }
  } else {
    GenericNodes *tnode;
    tnode = tnodes[0];
    for (auto i = 0; i < 3; ++i) {
      if (dofs(i) != 0) {
        tnode->setBoundaryCondition(i);
      }
    }
  }
}

void VertexRuntime::setLoad(LoadList &loadlist, indexType meshid,
                            ShapeFunctionTypes shapeType, indexType shapeOrder,
                            Types::VectorX<prec> &Loads, indexType propNumber,
                            Types::VectorX<prec> &direction, bool local,
                            bool add) {

  std::vector<GenericNodes *> nodeVector;
  this->getNodes(nodeVector, meshid);

  if (!nodeVector.empty()) {
    auto tempDofs = nodeVector[0]->getDegreesOfFreedom();

    indexType nn = Loads.size();

    for (auto j = 0; j < nn; ++j) {
      loadlist.setLoad(propNumber, tempDofs[j]->getId(),
                                      Loads(j), add);
    }
  }
}

void VertexRuntime::setPrescribedSolution(
    LoadList &loadlist, indexType meshid, ShapeFunctionTypes shapeType,
    indexType shapeOrder, Types::VectorX<prec> &Solution, indexType propNumber,
    Types::VectorX<prec> &direction, bool local, bool add) {
  std::vector<GenericNodes *> nodeVector;
  this->getNodes(nodeVector, meshid);

  if (!nodeVector.empty()) {
    auto tempDofs = nodeVector[0]->getDegreesOfFreedom();

    indexType nn = Solution.size();

    for (auto j = 0; j < nn; ++j) {
      if (Solution(j) != prec(0)) {
        loadlist.setLoad(
            propNumber, tempDofs[j]->getId(), Solution(j), add);
        // tempDofs[j]->setStatus(dofStatus::inactive);
      }
    }
  }
}

void VertexRuntime::getNodes(std::vector<GenericNodes *> &nodeVector,
                             indexType meshId) {
  nodeVector = this->getNodesOfSet(meshId);
}

void VertexRuntime::geometryToParaview(vtkPlotInterface &paraviewAdapter,
                                       indexType mainMesh, indexType subMesh) {
  Types::Vector3<prec> coor = m_Vertex_Data_Element->getCoordinates(); 
  paraviewAdapter.addPoint(mainMesh, subMesh, m_Data_Element->getId(),
                           coor);
}

void VertexRuntime::connectEdge(indexType edgeId) {
  m_Vertex_Data_Element->connectEdge(edgeId);
}

auto VertexRuntime::getConnectedEdges() -> std::vector<indexType> {
  return m_Vertex_Data_Element->getConnectedEdges();
}

void VertexRuntime::connectFace(indexType faceId) {
  m_Vertex_Data_Element->connectFace(faceId);
}

auto VertexRuntime::getConnectedFaces() -> std::vector<indexType> {
  return m_Vertex_Data_Element->getConnectedFaces();
}

const GeometryTypes VertexRuntime::type = GeometryTypes::Vertex;

} // namespace HierAMuS::Geometry
