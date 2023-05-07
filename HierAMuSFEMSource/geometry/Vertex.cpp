// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "equations/GenericNodes.h"
#include <equations/DegreeOfFreedom.h>

#include <geometry/GeometryTypes.h>
#include <geometry/Vertex.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <pointercollection/pointercollection.h>

#include <loads/LoadList.h>

#include <iomanip>
#include <vector>

#include <plot/vtkplotClass.h>

namespace HierAMuS::Geometry {

Vertex::Vertex()
    : connectedEdges({-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}),
      connectedFaces({-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}) {
  this->coors.setZero();
}

Vertex::~Vertex() {}

auto Vertex::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::Vertex::type;
}

void Vertex::setCoordinates(prec x, prec y, prec z) {
  this->coors(0) = x;
  this->coors(1) = y;
  this->coors(2) = z;
}

void Vertex::setCoordinates(const Types::Vector3<prec> &coorIn) {
  this->coors = coorIn;
}

inline void Vertex::print(PointerCollection &pointers) {
  pointers.getSPDLogger().debug(*this);
  this->printEqInfo(pointers);
}

//
// void Vertex::print(){
//	std::cout << "Vertex id: " << this->id  <<  " x: " << std::setw(8) <<
// this->x
//			<< " y: " << std::setw(8) << this->y << " z: " <<
// std::setw(8)<< this->z << std::endl;
//}

auto Vertex::getCoordinates() -> Types::Vector3<prec> { return this->coors; }

auto Vertex::getCoordinates(PointerCollection &pointers,
                            IntegrationPoint &IntPoint)
    -> Types::Vector3<prec> {
  return this->coors;
}

void Vertex::setBoundaryCondition(PointerCollection &pointers, indexType meshId,
                                  indexType order, ShapeFunctionTypes shapeType,
                                  Types::Vector3<indexType> &dofs, bool set) {
  std::vector<GenericNodes *> tnodes;
  this->getNodes(pointers, tnodes, meshId);

  if (set) {
    GenericNodes *tnode;
    tnode = tnodes[0];
    for (auto i = 0; i < 3; ++i) {
      if (dofs(i) != 0) {
        tnode->setBoundaryCondition(pointers, i);
      } else {
        tnode->unsetBoundaryCondition(pointers, i);
      }
    }
  } else {
    GenericNodes *tnode;
    tnode = tnodes[0];
    for (auto i = 0; i < 3; ++i) {
      if (dofs(i) != 0) {
        tnode->setBoundaryCondition(pointers, i);
      }
    }
  }
}

void Vertex::setLoad(PointerCollection &pointers, indexType meshid,
                     ShapeFunctionTypes shapeType, indexType shapeOrder,
                     Types::VectorX<prec> &Loads, indexType propNumber,
                     Types::VectorX<prec> &direction, bool local, bool add) {

  std::vector<GenericNodes *> nodeVector;
  this->getNodes(pointers, nodeVector, meshid);

  if (!nodeVector.empty()) {
    auto tempDofs = nodeVector[0]->getDegreesOfFreedom(pointers);

    indexType nn = Loads.size();

    for (auto j = 0; j < nn; ++j) {
      pointers.getLoadList()->setLoad(propNumber, tempDofs[j]->getId(),
                                      Loads(j), add);
    }
  }
}

void Vertex::setPrescribedSolution(
    PointerCollection &pointers, indexType meshid, ShapeFunctionTypes shapeType,
    indexType shapeOrder, Types::VectorX<prec> &Solution, indexType propNumber,
    Types::VectorX<prec> &direction, bool local, bool add) {
  std::vector<GenericNodes *> nodeVector;
  this->getNodes(pointers, nodeVector, meshid);

  if (!nodeVector.empty()) {
    auto tempDofs = nodeVector[0]->getDegreesOfFreedom(pointers);

    indexType nn = Solution.size();

    for (auto j = 0; j < nn; ++j) {
      if (Solution(j) != prec(0)) {
        pointers.getPrescribedDisplacements()->setLoad(
            propNumber, tempDofs[j]->getId(), Solution(j), add);
        //tempDofs[j]->setStatus(dofStatus::inactive);
      }
    }
  }
}

void Vertex::getNodes(PointerCollection &pointers,
                      std::vector<GenericNodes *> &nodeVector,
                      indexType meshId) {
  nodeVector.clear();
  this->getNodesOfSet(pointers, nodeVector, meshId);
}

void Vertex::geometryToParaview(PointerCollection &pointers,
                                vtkPlotInterface &paraviewAdapter,
                                indexType mainMesh, indexType subMesh) {
  paraviewAdapter.addPoint(mainMesh, subMesh, this->id, this->coors);
}

void Vertex::connectEdge(indexType edgeId) {

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

auto Vertex::getConnectedEdges() -> std::vector<indexType> {
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

void Vertex::connectFace(indexType faceId) {

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

auto Vertex::getConnectedFaces() -> std::vector<indexType> {
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

const GeometryTypes Vertex::type = GeometryTypes::Vertex;

} // namespace HierAMuS::Geometry
