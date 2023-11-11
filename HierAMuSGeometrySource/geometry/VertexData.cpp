// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "GenericNodes.h"

#include <geometry/GeometryTypes.h>
#include <geometry/VertexData.h>

#include "LoadList.h"

#include <vector>

#include <plot/vtkplotClass.h>

namespace HierAMuS::Geometry {

VertexData::VertexData()
    : connectedEdges({-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}),
      connectedFaces({-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}) {
  this->coors.setZero();
}

VertexData::~VertexData() {}

auto VertexData::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::VertexData::type;
}

void VertexData::setCoordinates(prec x, prec y, prec z) {
  this->coors(0) = x;
  this->coors(1) = y;
  this->coors(2) = z;
}

void VertexData::setCoordinates(const Types::Vector3<prec> &coorIn) {
  this->coors = coorIn;
}

void VertexData::print(spdlog::logger &Logger) {
  Logger.debug(*this);
  this->printEqInfo(Logger);
}


auto VertexData::getCoordinates() -> Types::Vector3<prec> { return this->coors; }



void VertexData::setBoundaryCondition(indexType meshId,
                                  indexType order, ShapeFunctionTypes shapeType,
                                  Types::Vector3<indexType> &dofs, bool set) {
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

void VertexData::setLoad(LoadList &loadlist, indexType meshid,
                     ShapeFunctionTypes shapeType, indexType shapeOrder,
                     Types::VectorX<prec> &Loads, indexType propNumber,
                     Types::VectorX<prec> &direction, bool local, bool add) {

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

void VertexData::setPrescribedSolution(
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
        //tempDofs[j]->setStatus(dofStatus::inactive);
      }
    }
  }
}

void VertexData::getNodes(std::vector<GenericNodes *> &nodeVector,
                      indexType meshId) {
  nodeVector = this->getNodesOfSet(meshId);
}

void VertexData::geometryToParaview(vtkPlotInterface &paraviewAdapter,
                                indexType mainMesh, indexType subMesh) {
  paraviewAdapter.addPoint(mainMesh, subMesh, this->id, this->coors);
}

void VertexData::connectEdge(indexType edgeId) {

  bool search = true;
  bool shift = false;
  indexType pos = 0;
  while (search) {
    if (this->connectedEdges[pos] == -1 ||
        this->connectedEdges[pos] == edgeId) {
      search = false;
    } else if (this->connectedEdges[pos] > edgeId) {
      shift = true;
      search = false;
    } else {
      pos++;
    }
  }
  if (shift) {
    search = true;
    indexType pos2 = pos;
    while (search) {
      if (this->connectedEdges[pos2] == -1) {
        search = false;
      } else {
        pos2++;
      }
    }
    for (indexType i = pos2; i > pos; i--) {
      this->connectedEdges[i] = this->connectedEdges[i - 1];
    }
    this->connectedEdges[pos] = edgeId;
  } else {
    this->connectedEdges[pos] = edgeId;
  }
}

auto VertexData::getConnectedEdges() -> std::vector<indexType> {
  bool search = true;
  indexType pos = 0;
  while (search) {
    if (this->connectedEdges[pos] == -1) {
      search = false;
    } else {
      pos++;
    }
  }
  // pos--;
  std::vector<indexType> ret;
  if (pos >= 0) {
    ret = std::vector<indexType>(this->connectedEdges.begin(),
                                 this->connectedEdges.begin() + pos);
  }

  return ret;
}

void VertexData::connectFace(indexType faceId) {

  bool search = true;
  bool shift = false;
  indexType pos = 0;
  while (search) {
    if (this->connectedFaces[pos] == -1 ||
        this->connectedFaces[pos] == faceId) {
      search = false;
    } else if (this->connectedFaces[pos] > faceId) {
      shift = true;
      search = false;
    } else {
      pos++;
    }
  }
  if (shift) {
    search = true;
    indexType pos2 = pos;
    while (search) {
      if (this->connectedFaces[pos2] == -1) {
        search = false;
      } else {
        pos2++;
      }
    }
    for (indexType i = pos2; i > pos; i--) {
      this->connectedFaces[i] = this->connectedFaces[i - 1];
    }
    this->connectedFaces[pos] = faceId;
  } else {
    this->connectedFaces[pos] = faceId;
  }
}

auto VertexData::getConnectedFaces() -> std::vector<indexType> {
  bool search = true;
  indexType pos = 0;
  while (search) {
    if (this->connectedFaces[pos] == -1) {
      search = false;
    } else {
      pos++;
    }
  }
  // pos--;
  std::vector<indexType> ret;
  if (pos >= 0) {
    ret = std::vector<indexType>(this->connectedFaces.begin(),
                                 this->connectedFaces.begin() + pos);
  }

  return ret;
}

const GeometryTypes VertexData::type = GeometryTypes::Vertex;

} // namespace HierAMuS::Geometry
