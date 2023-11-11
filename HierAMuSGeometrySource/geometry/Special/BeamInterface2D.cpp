// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include <types/MatrixTypes.h>

// Equations
#include "GenericNodes.h"


#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <geometry/Special/BeamInterface2D.h>
#include <pointercollection/pointercollection.h>

#include <geometry/Edges/EdgesData.h>
#include <geometry/GeometryData.h>
#include <geometry/VertexData.h>

#include <algorithm>
#include <iomanip>
#include <vector>

namespace HierAMuS::Geometry {

BeamInterface2D::BeamInterface2D() : Special(), initialized(false) {}

BeamInterface2D::~BeamInterface2D() = default;

auto BeamInterface2D::getGlobalEdgeNumber(indexType localEdgeNumber) -> indexType {
  return this->edges[localEdgeNumber];
}

auto BeamInterface2D::getLocalJacobian(PointerCollection &pointers,
                                       indexType localEdgeNumber, prec eta) -> prec {
  
  auto &edge = pointers.getGeometryData()->getEdgeData(this->edges[localEdgeNumber]);

  auto &V1 = *edge.getVertex(0);
  auto &V2 = *edge.getVertex(1);
  auto c2 = V2.getCoordinates();
  auto c1 = V1.getCoordinates();

  auto edir = c2 - c1;

  prec jacobi = edir.dot(this->A2) / prec(2);

  return jacobi;
}

auto BeamInterface2D::getCoordinates(PointerCollection &pointers,
                                     IntegrationPoint &IntPoint)
    -> Types::Vector3<prec>  {

  Types::Vector3<prec> coor;
  coor.setZero();
  auto &Edge =
      pointers.getGeometryData()->getEdgeData(this->edges[IntPoint.sectionNumber]);
  IntegrationPoint eInt;
  eInt.xi = IntPoint.eta;
  coor = Edge.getCoordinates(eInt);
  prec xip = IntPoint.xi + prec(1);
  prec dL = this->thickness / prec(2) * xip;
  coor += dL * this->A1;

  return coor;
}

void BeamInterface2D::setEdges(const std::vector<indexType> &edgesIn) {
  this->edges = edgesIn;
}

void BeamInterface2D::setBeamVertex(indexType beamVertIn) {
  this->beamNode = beamVertIn;
}

void BeamInterface2D::computeWarpingShapes(PointerCollection &pointers,
                                           indexType localOrder) {

  if (!this->initialized) {
    this->initialized = true;
    this->warpingOrder = localOrder;
    this->computeGeometry(pointers);
    this->setUpGlobalLocalVertexMapping(pointers);

    this->computeAlphaBetaParameters(pointers);
  }
}

void BeamInterface2D::getLocalWarpingShapesA1(
    PointerCollection &pointers, Types::VectorX<prec> &shapes,
    Types::VectorX<prec> &shapeDerivative, indexType localEdgeNumber,
    prec eta) {
  
  auto &ledge = pointers.getGeometryData()->getEdgeData(this->edges[localEdgeNumber]);

  Types::Vector3<prec> edgeDir = ledge.getDirectionVector();
  prec sign = edgeDir.dot(this->A2);
  sign < 0 ? sign = prec(-1) : sign = prec(1);

  prec jacobi = this->getLocalJacobian(pointers, localEdgeNumber, eta);
  Types::VectorX<prec> eshape;
  Types::VectorX<prec> eshapeDeriv;
  ledge.getH1Shapes(this->warpingOrder, eshape, eshapeDeriv, eta);
  eshapeDeriv /= jacobi;
  prec z = this->getCrossSectionPosition(pointers, eta, localEdgeNumber);

  shapes = this->alphaBeta.template topRows<1>();
  shapes += z * this->alphaBeta.template bottomRows<1>();
  shapeDerivative = this->alphaBeta.template bottomRows<1>();

  indexType tid = 0;
  indexType numVerts = ledge.getNumberOfVerts();
  indexType idd;
  for (auto i = 0; i < numVerts; ++i) {
    auto &Vert = *ledge.getVertex(i);
    idd = Vert.getId();
    idd = this->globalLocalVertIdMap[idd];
    shapes(idd) += eshape(i);
    shapeDerivative(idd) += eshapeDeriv(i);
    if (tid < idd)
      tid = idd;
  }
  idd = ledge.getId();
  idd = this->globalLocalEdgeIdMap[idd];
  for (auto i = 0; i < this->warpingOrder - 1; ++i) {
    shapes(idd) += eshape(i + numVerts);
    shapeDerivative(idd) += eshapeDeriv(i + numVerts);
    ++idd;
  }
}

void BeamInterface2D::getLocalWarpingShapesA2(
    PointerCollection &pointers, Types::VectorX<prec> &shapes,
    Types::VectorX<prec> &shapeDerivative, indexType localEdgeNumber,
    prec eta) {
  auto &ledge = pointers.getGeometryData()->getEdgeData(this->edges[localEdgeNumber]);

  Types::Vector3<prec> edgeDir = ledge.getDirectionVector();
  prec sign = edgeDir.dot(this->A2);
  sign < 0 ? sign = prec(-1) : sign = prec(1);

  prec jacobi = this->getLocalJacobian(pointers, localEdgeNumber, eta);
  Types::VectorX<prec> eshape;
  Types::VectorX<prec> eshapeDeriv;
  ledge.getH1Shapes(this->warpingOrder, eshape, eshapeDeriv, eta);
  eshapeDeriv /= jacobi;
  
  shapes = this->alphaA2;
  shapeDerivative.setZero();

  indexType tid = 0;
  indexType numVerts = ledge.getNumberOfVerts();
  indexType idd;
  for (auto i = 0; i < numVerts; ++i) {
    auto &Vert = *ledge.getVertex(i);
    idd = Vert.getId();
    idd = this->globalLocalVertIdMap[idd];
    shapes(idd) += eshape(i);
    shapeDerivative(idd) += eshapeDeriv(i);
    if (tid < idd)
      tid = idd;
  }
  idd = ledge.getId();
  idd = this->globalLocalEdgeIdMap[idd];
  for (auto i = 0; i < this->warpingOrder - 1; ++i) {
    shapes(idd) += eshape(i + numVerts);
    shapeDerivative(idd) += eshapeDeriv(i + numVerts);
    ++idd;
  }
}

void BeamInterface2D::getLocalH1ShapesA2(PointerCollection &pointers,
                                         Types::VectorX<prec> &shapes,
                                         Types::VectorX<prec> &shapeDerivative,
                                         indexType localEdgeNumber, prec eta,
                                         indexType shapeOrder,
                                         indexType meshId) {

  auto &edge = pointers.getGeometryData()->getEdgeData(this->edges[localEdgeNumber]);

  Types::VectorX<prec> eshape;
  Types::VectorX<prec> eshapeDeriv;
  edge.getH1Shapes(shapeOrder, eshape, eshapeDeriv, eta);

  prec dl = this->getLocalJacobian(pointers, localEdgeNumber, eta);
  eshapeDeriv /= dl;

  std::vector<GenericNodes *> tnodes;
  edge.getNodes(tnodes, meshId);
  shapes.resize(this->alphaBeta.cols());
  shapeDerivative.resize(this->alphaBeta.cols());
  shapes.setZero();
  shapeDerivative.setZero();

  std::map<indexType, indexType> nodemap;
  std::vector<GenericNodes *> surfaceNodes;
  this->getNodesOnSolid(pointers, surfaceNodes, meshId);
  nodemap = this->getNodeMapping(surfaceNodes);

  indexType numShapes;
  static_cast<indexType>(tnodes.size()) < shapeOrder + 1 ? numShapes = tnodes.size()
                                 : numShapes = shapeOrder + 1;

  indexType numVerts = edge.getNumberOfVerts();

  for (auto i = 0; i < numVerts; ++i) {
    auto &tVert = *edge.getVertex(i);
    indexType idd = tVert.getId();
    idd = this->globalLocalVertIdMap[idd];
    shapes(idd) = eshape(i);
    shapeDerivative(idd) = eshapeDeriv(i);
  }
  indexType idd = edge.getId();
  idd = this->globalLocalEdgeIdMap[idd];
  indexType cc = 0;
  for (auto i = 0; i < numShapes - numVerts; ++i) {
    shapes(idd) = eshape(i + numVerts + cc);
    shapeDerivative(idd) = eshapeDeriv(i + numVerts + cc);
    ++cc;
  }
}

void BeamInterface2D::getH1ShapesA1(PointerCollection &pointers,
                                    Types::VectorX<prec> &shapes,
                                    Types::VectorX<prec> &shapeDerivative,
                                    prec xi, indexType order) {

  auto &tedge = pointers.getGeometryData()->getEdgeData(this->edges[0]);
  tedge.getH1Shapes(order, shapes, shapeDerivative, xi);

  auto &Vert = pointers.getGeometryData()->getVertexData(this->beamNode);
  Types::Vector3<prec> ldir;
  ldir = Vert.getCoordinates() - this->surfaceCoordinat;
  prec dlen = ldir.dot(this->A1);
  dlen /= prec(2);
  shapeDerivative /= dlen;
}

void BeamInterface2D::globalToLocalEdgeNumbers(
    std::vector<indexType> &edgeNumbersInOut) {

  for (auto &i : edgeNumbersInOut) {
    auto itr = std::find(this->edges.begin(), this->edges.end(), i);
    if (itr == this->edges.end()) {
      throw std::runtime_error(
          "Error in beaminterface2D globalToLocalEdgeNumbers, not all edges in "
          "list");
}
    i = std::distance(this->edges.begin(), itr);
  }
}

void BeamInterface2D::computeGeometry(PointerCollection &pointers) {

  auto &edge = pointers.getGeometryData()->getEdgeData(this->edges[0]);
  auto &Vert = pointers.getGeometryData()->getVertexData(this->beamNode);

  Types::Vector3<prec> coorA;
  Types::Vector3<prec> coorB;
  Types::Vector3<prec> beamCoor;

  coorA = edge.getCoordinates(prec(-1));
  coorB = edge.getCoordinates(prec(1));
  beamCoor = Vert.getCoordinates();

  this->A2 = coorB - coorA;
  this->A2.normalize();

  this->A1(0) = this->A2(1);
  this->A1(1) = -this->A2(0);
  this->A1(2) = prec(0);

  Types::Vector3<prec> tempVec;
  tempVec = beamCoor - coorA;

  this->thickness = tempVec.dot(this->A1);

  if (this->thickness < prec(0)) {
    this->A1 *= prec(-1);
    this->thickness *= prec(-1);
    this->A2(0) = -this->A1(1);
    this->A2(1) = this->A1(0);
  }

  this->surfaceCoordinat = beamCoor - this->thickness * this->A1;
}

auto BeamInterface2D::getCrossSectionPosition(PointerCollection &pointers,
                                              prec eta, indexType edgenum)
    -> prec {
  auto &edge = pointers.getGeometryData()->getEdgeData(this->edges[edgenum]);
  Types::Vector3<prec> coor = edge.getCoordinates(eta);

  coor -= this->surfaceCoordinat;

  prec z = coor.dot(this->A2);

  return z;
}

auto BeamInterface2D::getDA(PointerCollection &pointers,
                            indexType localEdgeNumber, prec xi, prec eta) -> prec {

  prec dA;
  dA = this->thickness / prec(2);
  auto &tedge = pointers.getGeometryData()->getEdgeData(this->edges[localEdgeNumber]);
  auto &V1 = *tedge.getVertex(0);
  auto &V2 = *tedge.getVertex(1);
  auto c2 = V2.getCoordinates();
  auto c1 = V1.getCoordinates();
  auto dist = c2-c1;
  prec temp;
  temp = dist.norm();
  temp /= prec(2);
  dA *= temp;

  return dA;
}

void BeamInterface2D::setH1ShapesSurface(PointerCollection &pointers,
                                         indexType meshId, indexType order,
                                         NodeTypes type) {
  for (auto &i : this->edges) {
    auto &edge = pointers.getGeometryData()->getEdgeData(i);
    edge.setH1Shapes(meshId, order, type);
  }
}

void BeamInterface2D::setDofBeamNode(PointerCollection &pointers,
                                     indexType meshId, NodeTypes type) {
  auto &vert = pointers.getGeometryData()->getVertexData(this->beamNode);
  vert.setNodeSet(meshId, 1, type);
}

void BeamInterface2D::setUpGlobalLocalVertexMapping(
    PointerCollection &pointers) {

  
  indexType cc = 0;
  for (auto &i : this->edges) {
    auto &edge = pointers.getGeometryData()->getEdgeData(i);
    for (auto j = 0; j < edge.getNumberOfVerts(); ++j) {
      auto &vert = *edge.getVertex(j);
      indexType idd = vert.getId();
      if (this->globalLocalVertIdMap.find(idd) ==
          this->globalLocalVertIdMap.end()) {
        this->globalLocalVertIdMap[idd] = cc;
        ++cc;
      }
    }
  }
  for (auto &i : this->edges) {
    auto &edge = pointers.getGeometryData()->getEdgeData(i);
    indexType idd = edge.getId();
    if (this->globalLocalEdgeIdMap.find(idd) ==
        this->globalLocalEdgeIdMap.end()) {
      this->globalLocalEdgeIdMap[idd] = cc;
      cc += this->warpingOrder - 1;
    }
  }
  this->totalShapes = cc;
}

void BeamInterface2D::getDofsOnSolid(PointerCollection &pointers,
                                     std::vector<DegreeOfFreedom *> &Dofs,
                                     indexType meshID, indexType localOrder) {

  std::vector<DegreeOfFreedom *> tempDofs;
  
  Dofs.resize((this->edges.size() + 1) * 3);

  for (auto &i : this->edges) {
    auto &tedge = pointers.getGeometryData()->getEdgeData(i);
    indexType nVerts = tedge.getNumberOfVerts();
    for (auto j = 0; j < nVerts; ++j) {
      auto &Vert = *tedge.getVertex(j);
      indexType idd = Vert.getId();
      idd = this->globalLocalVertIdMap[idd];
      idd *= 3;
      Vert.getNodeEquationIds(tempDofs, meshID, 0);
      for (auto k = 0; k < tempDofs.size(); ++k) {
        Dofs[idd + k] = tempDofs[k];
      }
    }
  }
  //indexType startPos = this->globalLocalVertIdMap.size();
  //startPos *= 3;

  for (auto &i : this->edges) {
    auto &tedge = pointers.getGeometryData()->getEdgeData(i);
    tempDofs.clear();
    tedge.getH1DofsInternal(tempDofs, meshID, localOrder);
    Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());
  }

  // for (auto &i : this->globalLocalVertIdMap) {
  //  Vert = pointers.getGeometryData()->getVertex(i.first);
  //  Vert->getNodeEquationIds(pointers, tempDofs, meshID, 0);
  //  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());
  //}
  // for (auto& i : this->edges) {
  //    Edges* tedge;
  //    tedge = pointers.getGeometryData()->getEdge(i);
  //    tempDofs.clear();
  //    tedge->getH1DofsInternal(pointers, tempDofs, meshID, localOrder);
  //    Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());
  //}
}

void BeamInterface2D::getNodesOnSolid(PointerCollection &pointers,
                                      std::vector<GenericNodes *> &nodes,
                                      indexType meshId) {
  std::vector<GenericNodes *> tempNodes;

  std::map<indexType, indexType> tempMap;

  for (auto &i : this->edges) {
    auto &edge = pointers.getGeometryData()->getEdgeData(i);
    tempNodes.clear();
    edge.getNodes(tempNodes, meshId);
    for (auto &j : tempNodes) {
      indexType idd = j->getId();
      if (tempMap.find(idd) == tempMap.end()) {
        tempMap[idd] = 1;
        nodes.push_back(j);
      }
    }
  }
}

auto
BeamInterface2D::getNodeMapping(std::vector<GenericNodes *> &nodes) -> std::map<indexType, indexType> {
  std::map<indexType, indexType> mapping;

  indexType cc = 0;
  for (auto &i : nodes) {
    indexType idd = i->getId();
    if (mapping.find(idd) == mapping.end()) {
      mapping[idd] = cc;
      ++cc;
    }
  }

  return mapping;
}

void BeamInterface2D::getDofsOnBeamNode(PointerCollection &pointers,
                                        std::vector<DegreeOfFreedom *> &Dofs,
                                        indexType meshId) {
  auto &Vert = pointers.getGeometryData()->getVertexData(this->beamNode);
  Vert.getNodeEquationIds(Dofs, meshId, 0);
}

void BeamInterface2D::getAllDofs(PointerCollection &pointers,
                                 std::vector<DegreeOfFreedom *> &Dofs) {

  std::vector<DegreeOfFreedom *> tempDofs;
  for (auto &ee : this->edges) {
    auto &edge = pointers.getGeometryData()->getEdgeData(ee);
    for (auto j = 0; j < edge.getNumberOfVerts(); ++j) {
      auto &vert = *edge.getVertex(j);
      vert.getAllEquationsIds(tempDofs);
      Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());
      tempDofs.clear();
    }
    edge.getAllEquationsIds(tempDofs);
    Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());
    tempDofs.clear();
  }
  auto &vert = pointers.getGeometryData()->getVertexData(this->beamNode);
  vert.getAllEquationsIds(tempDofs);
  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());
}

void BeamInterface2D::computeAlphaBetaParameters(PointerCollection &pointers) {

  auto GP = this->getIntegrationPoints(-1);
  GP.setOrder(this->warpingOrder * 2);
  GP.setNumberOfSections(1);
  
  Types::VectorX<prec> shape;
  Types::VectorX<prec> shapeDerivative;
  Types::Matrix22<prec> lgs;
  lgs.setZero();

  this->alphaA2.resize(this->totalShapes);
  this->alphaBeta.resize(2, this->totalShapes);
  this->alphaBeta.setZero();

  for (auto edgeNum = 0; edgeNum < this->edges.size(); ++edgeNum) {
    auto &edge = pointers.getGeometryData()->getEdgeData(this->edges[edgeNum]);

    for (auto ngp = 0; ngp < GP.getTotalGP(); ++ngp) {
      prec zCoordinate =
          this->getCrossSectionPosition(pointers, GP.getXi(ngp), edgeNum);
      edge.getH1Shapes(this->warpingOrder, shape, shapeDerivative,
                        GP.getXi(ngp));
      prec jacobi = edge.getJacobian(GP.getXi(ngp));
      indexType id;
      indexType totVerts = edge.getNumberOfVerts();
      for (auto nvert = 0; nvert < totVerts; ++nvert) {
        id = edge.getVertex(nvert)->getId(); // get global id
        id = this->globalLocalVertIdMap[id];            // get local id
        this->alphaBeta(0, id) += shape(nvert) * jacobi * GP.getWeight(ngp);
        this->alphaBeta(1, id) +=
            zCoordinate * shape(nvert) * jacobi * GP.getWeight(ngp);
      }
      id = edge.getId();
      id = this->globalLocalEdgeIdMap[id];
      for (auto i = 0; i < this->warpingOrder - 1; ++i) {
        this->alphaBeta(0, id) +=
            shape(i + totVerts) * jacobi * GP.getWeight(ngp);
        this->alphaBeta(1, id) +=
            zCoordinate * shape(i + totVerts) * jacobi * GP.getWeight(ngp);
        ++id;
      }
      lgs(0, 0) += jacobi * GP.getWeight(ngp);
      lgs(1, 1) += zCoordinate * zCoordinate * jacobi * GP.getWeight(ngp);
      lgs(0, 1) += zCoordinate * jacobi * GP.getWeight(ngp);
    }
  }
  lgs(1, 0) = lgs(0, 1);
  this->alphaA2 = this->alphaBeta.template topRows<1>()/lgs(0,0);
  this->alphaBeta = -lgs.inverse() * this->alphaBeta;
  
}

void BeamInterface2D::setBCOnVert(PointerCollection &pointers,
                                  const indexType &meshID,
                                  const indexType &dof) {
  auto &Vert = pointers.getGeometryData()->getVertexData(this->beamNode);
  auto nodeList = Vert.getNodeSetNodeListMeshId(meshID);
  nodeList[0].setBoundaryCondition(dof);
}

void BeamInterface2D::setBCOnSolid(PointerCollection &pointers,
                                   indexType meshID, indexType dof) {
  std::vector<GenericNodes *> nodes;
  this->getNodesOnSolid(pointers, nodes, meshID);
  GenericNodes *tnode;
  for (auto &i : nodes) {
    tnode = i;
    tnode->setBoundaryCondition(dof);
  }
}

auto BeamInterface2D::getVertex(PointerCollection &pointers) -> VertexData & {
  return pointers.getGeometryData()->getVertexData(this->beamNode);
}

auto BeamInterface2D::getIntegrationPoints(indexType elementId) -> IntegrationPoints {
  auto GP = IntegrationPointsManagement::getIntegrationsPoints(this->id);
  GP.setType(Gauss2D);
  GP.setNumberOfSections(this->edges.size());
  return GP;
}

const GeometryTypes BeamInterface2D::type = GeometryTypes::BeamInterface2D;

} /* namespace HierAMuS::Geometry */

