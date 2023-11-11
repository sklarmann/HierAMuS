// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <iostream>
#include <map>


#include "beamInterfaceElement2D.h"
#include "datatypes.h"

#include "plot/vtkplotClassBase.h"
#include <types/MatrixTypes.h>

#include <finiteElements/GenericFiniteElement.h>
#include <finiteElements/beamInterfaceElement2D.h>

#include <materials/Material.h>

#include <elementFormulations/GenericElementFormulation.h>

#include <pointercollection/pointercollection.h>

#include <geometry/GeometryBaseData.h>
#include <geometry/Special/BeamInterface2D.h>
#include <geometry/Edges/EdgesData.h>
#include <geometry/VertexData.h>


#include <geometry/GeometryData.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>

#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "shapefunctions/IntegrationsPoints/dataClasses/GaussPoints.h"

#include <shapefunctions/LagrangeShape.h>
#include <solver/GenericSolutionState.h>
#include <solver/SolutionTypes.h>

#include <algorithm>
#include <sys/types.h>

#include "GenericNodes.h"



namespace HierAMuS::FiniteElement {
  

beamInterfaceElement2D::~beamInterfaceElement2D() = default;

void beamInterfaceElement2D::setVerts(std::vector<indexType> &vertsIn) {
  if (vertsIn.size() == 1) {
    this->vert = vertsIn[0];
  } else {
  }
}

auto
beamInterfaceElement2D::getGlobalEdgeNumber(PointerCollection& pointers, indexType localEdgeNumber) -> indexType {

  auto elem = this->getInterfaceGeoElem(pointers);
  return elem->getGlobalEdgeNumber(localEdgeNumber);
}

void beamInterfaceElement2D::getShapes(
    ptrCol &pointers, const indexType &order, const indexType &meshIDDisp,
    const indexType &meshIDWarp, const indexType &meshIDRot,
    Types::VectorX<prec> &shapeX, Types::Matrix2X<prec> &dshapeX,
    Types::VectorX<prec> &shapeY, Types::Matrix2X<prec> &dshapeY, prec &detj,
    const prec &xi, const prec &eta, const indexType &localedgeNumber) {
  auto &tedge =
      pointers.getGeometryData()->getEdgeData(this->edges[localedgeNumber]);
  auto tvert = pointers.getGeometryData()->getVertexData(this->vert);
  Types::VectorX<prec> edgeShape;
  Types::VectorX<prec> dedgeShape;
  tedge.getH1Shapes(1, edgeShape, dedgeShape, eta);
  Types::Vector3<prec> ecoor1;
  Types::Vector3<prec> ecoor2;
  Types::Vector3<prec> tempcoor;
  ecoor1 = tedge.getCoordinates(prec(-1));
  ecoor2 = tedge.getCoordinates(prec(1));
  prec edetj = (ecoor2 - ecoor1).norm() / prec(2);
  prec sdetj = this->thickness / prec(2);
  ecoor1 = tedge.getCoordinates(eta);
  ecoor1 += this->thickness * this->normal;
  tempcoor = ecoor1 - tvert.getCoordinates();
  prec curry = tempcoor.dot(this->tangent);

  indexType nnodes = this->dofmapWarp.size();
  indexType totnodes = nnodes * 2 + 2;
  shapeX.resize(totnodes);
  dshapeX.resize(2, totnodes);
  shapeX.setZero();
  dshapeX.setZero();

  shapeY.resize(totnodes);
  dshapeY.resize(2, totnodes);
  shapeY.setZero();
  dshapeY.setZero();

  Types::VectorX<prec> tshp;
  Types::VectorX<prec> dtshp;
  Types::VectorX<prec> tshpx;
  Types::VectorX<prec> dtshpx;
  tedge.getH1Shapes(order, tshp, dtshp, eta);
  tedge.getH1Shapes(order, tshpx, dtshpx, xi);
  dtshpx /= sdetj;
  dtshp /= edetj;
  for (auto i = 0; i < nnodes; ++i) {
    shapeX(nnodes + i) =
        (this->paramsWarpingX(0, i) + curry * this->paramsWarpingX(1, i)) *
        tshpx(1);
    dshapeX(0, nnodes + i) =
        (this->paramsWarpingX(0, i) + curry * this->paramsWarpingX(1, i)) *
        dtshpx(1);
    // dshapeX(1, nnodes + i) = (edetj * this->paramsWx(1, i)) * tshpx(1);
    dshapeX(1, nnodes + i) = (this->paramsWarpingX(1, i)) * tshpx(1);

    shapeY(nnodes + i) = (this->paramsWarpingY(0, i)) * tshpx(1);
    // dshapeY(0, nnodes + i) = (this->paramsWy(0, i))* dtshpx(1) * tshpx(1);
    dshapeY(1, nnodes + i) =
        (this->paramsWarpingY(0, i)); //* dtshpx(1);// * tshpx(1);
    // if (curry < prec(0)) dshapeX(1, nnodes + i) *= prec(-1);
  }
  std::vector<GenericNodes *> tempNodes;
  tedge.getNodes(tempNodes, meshIDDisp);
  for (auto i = 0; i < tempNodes.size(); ++i) {
    indexType id = this->dofmapDisp[tempNodes[i]->getId()];
    shapeX(id) += tshp(i) * tshpx(0);
    dshapeX(0, id) += tshp(i) * dtshpx(0);
    dshapeX(1, id) += dtshp(i) * tshpx(0);

    shapeY(id) += tshp(i) * tshpx(0);
    dshapeY(0, id) += tshp(i) * dtshpx(0);
    dshapeY(1, id) += dtshp(i) * tshpx(0);
  }
  tedge.getNodes(tempNodes, meshIDWarp);
  for (auto i = 0; i < tempNodes.size(); ++i) {
    indexType id = this->dofmapWarp[tempNodes[i]->getId()] + nnodes;
    shapeX(id) += tshp(i) * tshpx(1);
    dshapeX(0, id) += tshp(i) * dtshpx(1);
    dshapeX(1, id) += dtshp(i) * tshpx(1);

    shapeY(id) += tshp(i) * tshpx(1);
    dshapeY(0, id) += tshp(i) * dtshpx(1);
    dshapeY(1, id) += tshp(i) * tshpx(1);
  }
  indexType bid = nnodes * indexType(2);
  shapeX(bid) = tshpx(1);
  shapeY(bid) = tshpx(1);
  dshapeX(0, bid) = dtshpx(1);
  dshapeY(0, bid) = dtshpx(1);

  shapeX(bid + 1) = -tshpx(1) * curry;
  dshapeX(0, bid + 1) = -dtshpx(1) * curry;
  // dshapeX(1, bid + 1) = -tshpx(1) * edetj;
  dshapeX(1, bid + 1) = -tshpx(1);

  Types::Matrix22<prec> jac;
  jac.setZero();
  jac(0, 0) = sdetj;
  jac(1, 1) = edetj;
  // dshapeX = (jac.inverse()).transpose() * dshapeX;
  // dshapeY = (jac.inverse()).transpose() * dshapeY;

  detj = edetj * sdetj;
}

void beamInterfaceElement2D::getEtaShape(
    ptrCol &pointers, const indexType &meshID, const indexType &edgeNum,
    const indexType &order, const prec &eta, Types::VectorX<prec> &shape,
    Types::VectorX<prec> &dshape) {
  indexType numNodes = this->dofmapDisp.size();
  shape.resize(numNodes);
  dshape.resize(numNodes);
  shape.setZero();
  dshape.setZero();

  auto &tedge = pointers.getGeometryData()->getEdgeData(this->edges[edgeNum]);
  Types::VectorX<prec> shp;
  Types::VectorX<prec> dshp;
  tedge.getH1Shapes(order, shp, dshp, eta);
  Types::Vector3<prec> coor1;
  Types::Vector3<prec> coor2;
  coor1 = tedge.getCoordinates(prec(-1));
  coor2 = tedge.getCoordinates(prec(1));
  coor2 -= coor1;

  prec len = coor2.norm();
  dshp *= prec(2) / len;

  std::vector<GenericNodes *> tempNodes;
  tedge.getNodes(tempNodes, meshID);

  for (auto i = 0; i < tempNodes.size(); ++i) {
    indexType id = this->dofmapDisp[tempNodes[i]->getId()];
    shape(id) = shp(i);
    dshape(id) = dshp(i);
  }
}

auto beamInterfaceElement2D::getDA(ptrCol &pointers, const indexType &edgeNum,
                                   const prec &xi, const prec &eta) -> prec {

  Geometry::BeamInterface2D *elem;
  elem = this->getInterfaceGeoElem(pointers);

  return elem->getDA(pointers, edgeNum, xi, eta);

  // Geometry::Edges *tedge =
  //     pointers.getGeometryData()->getEdge(this->edges[edgeNum]);
  // Types::Vector3<prec> coor1;
  // Types::Vector3<prec> coor2;
  // coor1 = tedge->getCoordinates(pointers, prec(-1));
  // coor2 = tedge->getCoordinates(pointers, prec(1));
  // coor2 -= coor1;

  // prec len = coor2.norm();

  // prec dA = len * this->thickness / prec(4);
  // return dA;
}

void beamInterfaceElement2D::computeShapes(ptrCol &pointers,
                                           const indexType &mesIDDisp,
                                           const indexType &mesIDWarp,
                                           const indexType &mesIDRot,
                                           const indexType &order) {
  // Constrain Dofs
  pointers.getGeometryData()
      ->getVertexData(this->vert)
      .setAllNodeBoundaryConditionMeshId(mesIDDisp, 2);
  pointers.getGeometryData()
      ->getVertexData(this->vert)
      .setAllNodeBoundaryConditionMeshId(mesIDRot, 1);
  pointers.getGeometryData()
      ->getVertexData(this->vert)
      .setAllNodeBoundaryConditionMeshId(mesIDRot, 2);

  std::vector<Geometry::VertexData *> tempElems;

  for (auto &i : this->edges) {
    auto &tempEdge = pointers.getGeometryData()->getEdgeData(i);
    tempEdge.getVerts(tempElems);
    tempEdge.setAllNodeBoundaryConditionMeshId(mesIDDisp, 2);
    tempEdge.setAllNodeBoundaryConditionMeshId(mesIDWarp, 2);
  }
  {
    auto &tempEdge = pointers.getGeometryData()->getEdgeData(this->edges[0]);
    tempEdge.getVerts(tempElems);
    tempElems[0]->setAllNodeBoundaryConditionMeshId(mesIDWarp, 0);
  }
  // tempElems[0]->setAllNodeBoundaryConditionMeshId(pointers, mesIDWarp, 1);
  // tempElems[1]->setAllNodeBoundaryConditionMeshId(pointers, mesIDWarp, 2);
  {
    auto &tempEdge = pointers.getGeometryData()->getEdgeData(
        this->edges[this->edges.size() - 1]);
    tempEdge.getVerts(tempElems);
    tempElems[1]->setAllNodeBoundaryConditionMeshId(mesIDWarp, 0);
  }

  std::vector<GenericNodes *> nodes;
  GenericNodes *tempNode;
  indexType nodeID;
  indexType nnodes = 0;
  indexType nnodes2 = 0;
  // ID Mapping of Warping nodes
  for (auto &i : this->edges) {
    auto &tempEdge = pointers.getGeometryData()->getEdgeData(i);
    tempEdge.getNodes(nodes, mesIDWarp);
    for (auto &j : nodes) {
      tempNode = j;
      nodeID = tempNode->getId();
      if (this->dofmapWarp.find(nodeID) == this->dofmapWarp.end()) {
        this->dofmapWarp[nodeID] = nnodes;
        ++nnodes;
      }
    }
    tempEdge.getNodes(nodes, mesIDDisp);
    for (auto &j : nodes) {
      tempNode = j;
      nodeID = tempNode->getId();
      if (this->dofmapDisp.find(nodeID) == this->dofmapDisp.end()) {
        this->dofmapDisp[nodeID] = nnodes2;
        ++nnodes2;
      }
    }
  }

  // Computation of warping shape functions
  Types::VectorX<prec> shape;
  Types::VectorX<prec> shapeDerivative;
  Types::Matrix22<prec> sys1;
  Types::Matrix33<prec> trafo;
  trafo.setZero();
  trafo.block(0, 0, 3, 1) = this->normal;
  trafo.block(0, 1, 3, 1) = this->tangent;
  // std::cout << trafo << std::endl;
  Types::Matrix2X<prec> rhs1;
  Types::Matrix2X<prec> rhs2;
  sys1.setZero();
  rhs1.resize(2, this->dofmapWarp.size());
  rhs1.setZero();
  std::vector<prec> gp;
  std::vector<prec> weight;
  linearGP(gp, weight, 3);
  prec jac;
  Types::Vector3<prec> ycoor;
  Types::Vector3<prec> refcoor;
  refcoor = pointers.getGeometryData()->getVertexData(this->vert).getCoordinates();
  for (auto &edge : this->edges) {
    auto &tempEdge = pointers.getGeometryData()->getEdgeData(edge);
    tempEdge.getNodes(nodes, mesIDWarp);
    for (auto i = 0; i < gp.size(); ++i) {
      tempEdge.getH1Shapes(order, shape, shapeDerivative, gp[i]);
      this->getJacobianInterface(pointers, jac, gp[i], edge);
      ycoor = tempEdge.getCoordinates(gp[i]);
      ycoor += this->thickness * this->normal;
      ycoor -= refcoor;
      ycoor = trafo * ycoor;

      sys1(0, 0) += jac * weight[i];
      sys1(0, 1) += ycoor(1) * jac * weight[i];
      sys1(1, 0) += ycoor(1) * jac * weight[i];
      sys1(1, 1) += ycoor(1) * ycoor(1) * jac * weight[i];
      indexType dofnum = 0;
      for (auto &nn : nodes) {
        indexType id = this->dofmapWarp[nn->getId()];
        rhs1(0, id) += jac * weight[i] * shape(dofnum);
        rhs1(1, id) += jac * weight[i] * shape(dofnum) * ycoor(1);

        ++dofnum;
      }
    }
  }
  this->paramsWarpingY =
      -rhs1.block(0, 0, 1, this->dofmapWarp.size()) / sys1(0, 0);
  this->paramsWarpingX = -sys1.inverse() * rhs1;
}

void beamInterfaceElement2D::computeShapes2(ptrCol &pointers,
                                            const indexType &mesIDDisp,
                                            const indexType &mesIDRot,
                                            const indexType &order) {
  // Constrain Dofs
  pointers.getGeometryData()
      ->getVertexData(this->vert)
      .setAllNodeBoundaryConditionMeshId(mesIDDisp, 2);
  pointers.getGeometryData()
      ->getVertexData(this->vert)
      .setAllNodeBoundaryConditionMeshId(mesIDRot, 1);
  pointers.getGeometryData()
      ->getVertexData(this->vert)
      .setAllNodeBoundaryConditionMeshId(mesIDRot, 2);

  std::vector<Geometry::VertexData *> tempElems;
  for (auto &i : this->edges) {
    auto &tempEdge = pointers.getGeometryData()->getEdgeData(i);
    tempEdge.getVerts(tempElems);
    tempEdge.setAllNodeBoundaryConditionMeshId(mesIDDisp, 2);
    for (auto &j : tempElems) {
      j->setAllNodeBoundaryConditionMeshId(mesIDDisp, 2);
    }
  }

  std::vector<GenericNodes *> nodes;
  GenericNodes *tempNode;
  indexType nodeID;

  indexType nnodes2 = 0;
  // ID Mapping of Warping nodes
  for (auto &i : this->edges) {
    auto &tempEdge = pointers.getGeometryData()->getEdgeData(i);
    tempEdge.getNodes(nodes, mesIDDisp);
    for (auto &j : nodes) {
      tempNode = j;
      nodeID = tempNode->getId();
      if (this->dofmapDisp.find(nodeID) == this->dofmapDisp.end()) {
        this->dofmapDisp[nodeID] = nnodes2;
        ++nnodes2;
      }
    }
  }

  // Computation of warping shape functions
  Types::VectorX<prec> shape;
  Types::VectorX<prec> shapeDerivative;
  Types::Matrix22<prec> sys1;

  Types::Matrix33<prec> trafo;
  trafo.setZero();
  trafo.block(0, 0, 3, 1) = this->normal;
  trafo.block(0, 1, 3, 1) = this->tangent;

  Types::Matrix2X<prec> rhs1;
  Types::Matrix2X<prec> rhs2;
  sys1.setZero();
  rhs1.resize(2, this->dofmapDisp.size());
  rhs1.setZero();
  std::vector<prec> gp;
  std::vector<prec> weight;
  linearGP(gp, weight, 2);
  prec jac;
  Types::Vector3<prec> ycoor;
  Types::Vector3<prec> refcoor;
  refcoor = pointers.getGeometryData()->getVertexData(this->vert).getCoordinates();
  for (auto &edge : this->edges) {
    auto &tempEdge = pointers.getGeometryData()->getEdgeData(edge);
    tempEdge.getNodes(nodes, mesIDDisp);
    for (auto i = 0; i < gp.size(); ++i) {
      tempEdge.getH1Shapes(order, shape, shapeDerivative, gp[i]);
      this->getJacobianInterface(pointers, jac, gp[i], edge);
      ycoor = tempEdge.getCoordinates(gp[i]);
      ycoor += this->thickness * this->normal;
      ycoor -= refcoor;
      ycoor = trafo * ycoor;

      sys1(0, 0) += jac * weight[i];
      // sys1(0, 1) += ycoor(1) * jac * weight[i];
      // sys1(1, 0) += ycoor(1) * jac * weight[i];
      sys1(1, 1) += ycoor(1) * ycoor(1) * jac * weight[i];
      indexType dofnum = 0;
      for (auto &nn : nodes) {
        indexType id = this->dofmapDisp[nn->getId()];
        rhs1(0, id) += jac * weight[i] * shape(dofnum);
        rhs1(1, id) += jac * weight[i] * shape(dofnum) * ycoor(1);

        ++dofnum;
      }
    }
  }
  this->paramsWarpingX = sys1.inverse() * rhs1;
}

void beamInterfaceElement2D::setH1Shapes(ptrCol &pointers, indexType meshID,
                                         indexType order) {}

void beamInterfaceElement2D::getXiShape(ptrCol &pointers,
                                        const indexType &order, const prec &xi,
                                        Types::VectorX<prec> &shape,
                                        Types::VectorX<prec> &dshape) {

  Geometry::BeamInterface2D *elem;
  elem = this->getInterfaceGeoElem(pointers);
  elem->getH1ShapesA1(pointers, shape, dshape, xi, order);

  // Geometry::Edges *tedge;
  // tedge = pointers.getGeometryData()->getEdge(this->edges[0]);

  // tedge->getH1Shapes(pointers, order, shape, dshape, xi);
  // dshape *= prec(2);
  // dshape /= this->thickness;
}

auto beamInterfaceElement2D::getZCoordinate(ptrCol &pointers,
                                            const indexType &edgeNum,
                                            const prec &eta) -> prec {

  Geometry::BeamInterface2D *elem;
  elem = this->getInterfaceGeoElem(pointers);
  return elem->getCrossSectionPosition(pointers, eta, edgeNum);
}

void beamInterfaceElement2D::setNodeMapDisp(ptrCol &pointers,
                                            const indexType &meshId) {
  std::vector<GenericNodes *> Nodes;
  this->getNodesOnSolid(pointers, Nodes, meshId);
  GenericNodes *tnode;
  indexType pos = 0;
  for (auto & Node : Nodes) {
    tnode = Node;
    indexType nid = tnode->getId();
    if (this->NodeMapDisp.find(nid) == this->NodeMapDisp.end()) {
      this->NodeMapDisp[nid] = pos;
      ++pos;
    }
  }
}

void beamInterfaceElement2D::setNodeMapWarp(ptrCol &pointers,
                                            const indexType &meshId) {
  std::vector<GenericNodes *> Nodes;
  this->getNodesOnSolid(pointers, Nodes, meshId);
  GenericNodes *tnode;
  indexType pos = 0;
  for (auto & Node : Nodes) {
    tnode = Node;
    indexType nid = tnode->getId();
    if (this->NodeMapWarp.find(nid) == this->NodeMapWarp.end()) {
      this->NodeMapWarp[nid] = pos;
      ++pos;
    }
  }
}

void beamInterfaceElement2D::setDofsOnSolid(ptrCol &pointers,
                                            const indexType &meshID,
                                            const indexType &order) {

  Geometry::BeamInterface2D *ielem;

  ielem = this->getInterfaceGeoElem(pointers);

  ielem->setH1ShapesSurface(pointers, meshID, order, NodeTypes::displacement);
}

void beamInterfaceElement2D::setDofsOnVert(ptrCol &pointers,
                                           const indexType &meshID) {

  Geometry::BeamInterface2D *ielem;

  ielem = this->getInterfaceGeoElem(pointers);
  ielem->setDofBeamNode(pointers, meshID, NodeTypes::displacement);
}

void beamInterfaceElement2D::getDofsOnSolid(
    ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
    const indexType &meshID, indexType order) {

  Geometry::BeamInterface2D *elem;
  elem = this->getInterfaceGeoElem(pointers);
  elem->getDofsOnSolid(pointers, Dofs, meshID, order);

  // Dofs.clear();
  // std::vector<GenericNodes *> orderedNodeList;

  // this->getNodesOnSolid(pointers, orderedNodeList, meshID);
  // std::vector<DegreeOfFreedom *> tdofs;
  // for (auto &i : orderedNodeList) {
  //   i->getDegreesOfFreedom(pointers, tdofs);
  //   Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  // }
}

void beamInterfaceElement2D::getNodesOnSolid(ptrCol &pointers,
                                             std::vector<GenericNodes *> &Nodes,
                                             const indexType &meshID) {
  Nodes.clear();
  Geometry::BeamInterface2D *elem;
  elem = this->getInterfaceGeoElem(pointers);
  elem->getNodesOnSolid(pointers, Nodes, meshID);
}

void beamInterfaceElement2D::computeShapesLocalWarping(ptrCol &pointers,
                                                       const indexType &meshid,
                                                       const indexType &order) {

  // Get Integration Points, set type and order
  IntegrationPoints GP = this->getIntegrationPoints(pointers);
  GP.setOrder(order);
  GP.setNumberOfSections(1);

  // Computation of warping shape functions
  Types::VectorX<prec> shape;
  Types::VectorX<prec> shapeDerivative;
  Types::Matrix22<prec> sys1;
  Types::Matrix33<prec> trafo;
  trafo.setZero();
  trafo.block(0, 0, 3, 1) = this->normal;
  trafo.block(0, 1, 3, 1) = this->tangent;

  Types::Matrix2X<prec> rhs1;
  Types::VectorX<prec> rhs2;
  sys1.setZero();
  rhs1.resize(2, this->NodeMapWarp.size());
  rhs1.setZero();
  rhs2.resize(this->NodeMapWarp.size());
  rhs2.setZero();

  prec jac;
  Types::Vector3<prec> ycoor;
  Types::Vector3<prec> refcoor;
  refcoor = pointers.getGeometryData()->getVertexData(this->vert).getCoordinates();
  std::vector<GenericNodes *> nodes;
  indexType edgeCounter = 0;
  for (auto &edge : this->edges) {
    auto &tempEdge = pointers.getGeometryData()->getEdgeData(edge);
    tempEdge.getNodes(nodes, meshid);
    for (auto i = 0; i < GP.getTotalGP(); ++i) {
      tempEdge.getH1Shapes(order, shape, shapeDerivative,
                           GP.getXi(i));
      // jac = this->getJacEta(pointers, edge);
      this->getJacobianInterface(pointers, jac, GP.getXi(i), edge);

      prec y = this->getZCoordinate(pointers, edgeCounter, GP.getXi(i));

      sys1(0, 0) += jac * GP.getWeight(i);
      sys1(0, 1) += y * jac * GP.getWeight(i);
      sys1(1, 0) += y * jac * GP.getWeight(i);
      sys1(1, 1) += y * y * jac * GP.getWeight(i);
      indexType dofnum = 0;
      for (auto &nn : nodes) {
        indexType id = this->NodeMapWarp[nn->getId()];
        rhs1(0, id) += jac * GP.getWeight(i) * shape(dofnum);
        rhs1(1, id) += jac * GP.getWeight(i) * shape(dofnum) * y;

        rhs2(id) += jac * GP.getWeight(i) * shape(dofnum);

        ++dofnum;
      }
    }
    ++edgeCounter;
  }
  this->paramsWarpingX = -sys1.inverse() * rhs1;
  this->paramsWarpingY = -rhs2 / sys1(0, 0);
}

void beamInterfaceElement2D::computeShapesPolynomialWarping(
    ptrCol &pointers, const indexType &meshid, const indexType &order) {}

void beamInterfaceElement2D::getWarpingShapes(
    ptrCol &pointers, Types::VectorX<prec> &shapeWx,
    Types::VectorX<prec> &shapeWy, Types::VectorX<prec> &shapeWxy,
    Types::VectorX<prec> &shapeWyy, indexType edgeNum, prec eta,
    indexType order, indexType meshId) {

  shapeWx.resize(this->NodeMapWarp.size());
  shapeWy.resize(this->NodeMapWarp.size());
  shapeWxy.resize(this->NodeMapWarp.size());
  shapeWyy.resize(this->NodeMapWarp.size());

  shapeWx.setZero();
  shapeWy.setZero();
  shapeWxy.setZero();
  shapeWyy.setZero();

  auto &tEdge = pointers.getGeometryData()->getEdgeData(this->edges[edgeNum]);

  prec y = this->getZCoordinate(pointers, edgeNum, eta);

  shapeWx =
      this->paramsWarpingX.block(0, 0, 1, this->NodeMapWarp.size()).transpose();
  shapeWx +=
      y *
      this->paramsWarpingX.block(1, 0, 1, this->NodeMapWarp.size()).transpose();
  shapeWxy =
      this->paramsWarpingX.block(1, 0, 1, this->NodeMapWarp.size()).transpose();

  shapeWy = this->paramsWarpingY;

  Types::Vector3<prec> c1;
  Types::Vector3<prec> c2;
  c1 = tEdge.getCoordinates(-prec(1));
  c2 = tEdge.getCoordinates(prec(1));
  auto dir = c2 - c1;

  prec jacEdge = dir.transpose() * this->tangent;
  jacEdge /= prec(2);

  Types::VectorX<prec> shapeEdge;
  Types::VectorX<prec> diffShapeEdge;
  tEdge.getH1Shapes(order, shapeEdge, diffShapeEdge, eta);
  diffShapeEdge /= jacEdge;

  std::vector<GenericNodes *> nodes;
  tEdge.getNodes(nodes, meshId);
  indexType cc = 0;
  for (auto &nn : nodes) {
    indexType pos = this->NodeMapWarp[nn->getId()];
    shapeWxy(pos) += diffShapeEdge(cc);
    shapeWyy(pos) += diffShapeEdge(cc);

    shapeWx(pos) += shapeEdge(cc);
    shapeWy(pos) += shapeEdge(cc);
    ++cc;
  }
  // std::cout << shapeWx << std::endl;
}

void beamInterfaceElement2D::computeSurfaceDispShapes(ptrCol &pointers,
                                                      const indexType &meshID,
                                                      const indexType &order) {

  // Get Integration Points, set type and order
  IntegrationPoints GP = this->getIntegrationPoints(pointers);
  GP.setOrder(order + 1);
  GP.setNumberOfSections(1);

  // Computation of warping shape functions
  Types::VectorX<prec> shape;
  Types::VectorX<prec> shapeDerivative;
  Types::Matrix33<prec> trafo;
  trafo.setZero();
  trafo.block(0, 0, 3, 1) = this->normal;
  trafo.block(0, 1, 3, 1) = this->tangent;

  this->Nu.resize(this->NodeMapDisp.size());
  this->Nbeta.resize(this->NodeMapDisp.size());
  this->Nu.setZero();
  this->Nbeta.setZero();

  auto geoData = pointers.getGeometryData();
  std::vector<GenericNodes *> nodes;

  prec jac;
  prec A = prec(0);
  prec Iz = prec(0);
  Types::Vector3<prec> ycoor;
  for (auto &edge : this->edges) {
    auto &tedge = geoData->getEdgeData(edge);
    tedge.getNodes(nodes, meshID);
    for (auto i = 0; i < GP.getTotalGP(); ++i) {
      tedge.getH1Shapes(order, shape, shapeDerivative, GP.getXi(i));
      this->getJacobianInterface(pointers, jac, GP.getXi(i), edge);
      ycoor = tedge.getCoordinates(GP.getXi(i));
      ycoor -= this->refCoordinates;
      prec y = ycoor.transpose() * this->tangent;

      A += jac * GP.getWeight(i);
      Iz += y * y * jac * GP.getWeight(i);
      indexType dofnum = 0;
      for (auto &nn : nodes) {
        indexType id = this->NodeMapDisp[nn->getId()];
        this->Nu(id) += jac * GP.getWeight(i) * shape(dofnum);
        this->Nbeta(id) += jac * GP.getWeight(i) * shape(dofnum) * y;
        ++dofnum;
      }
    }
  }
  this->Nu /= A;
  this->Nbeta /= Iz;
}

void beamInterfaceElement2D::getSurfaceDispShapes(
    Types::VectorX<prec> &shapeNu, Types::VectorX<prec> &shapeNbeta) {
  shapeNu = this->Nu;
  shapeNbeta = this->Nbeta;
}

void beamInterfaceElement2D::getLocalSurfaceDispShapesSorted(
    ptrCol &pointers, Types::VectorX<prec> &shape, Types::VectorX<prec> &dshape,
    indexType edgeNum, prec eta, indexType order, indexType meshId) {

  Geometry::BeamInterface2D *elem;
  elem = this->getInterfaceGeoElem(pointers);
  elem->getLocalH1ShapesA2(pointers, shape, dshape, edgeNum, eta, order,
                           meshId);

  // Geometry::Edges *tEdge =
  //     pointers.getGeometryData()->getEdge(this->edges[edgeNum]);

  // auto c1 = tEdge->getCoordinates(pointers, -prec(1));
  // auto c2 = tEdge->getCoordinates(pointers, prec(1));

  // auto diff = c2 - c1;
  // prec len = diff.transpose() * this->tangent;
  // len /= prec(2);

  // Types::VectorX<prec> shapeEdge;
  // Types::VectorX<prec> dshapeEdge;
  // tEdge->getH1Shapes(pointers, order, shapeEdge, dshapeEdge, eta);
  // dshapeEdge /= len;

  // std::vector<GenericNodes *> Nodes;
  // tEdge->getNodes(pointers, Nodes, meshId);

  // shape.resize(this->NodeMapDisp.size());
  // dshape.resize(this->NodeMapDisp.size());
  // shape.setZero();
  // dshape.setZero();

  // indexType pos = 0;
  // for (auto &i : Nodes) {
  //   indexType id = this->NodeMapDisp[i->getId()];
  //   shape(id) = shapeEdge(pos);
  //   dshape(id) = dshapeEdge(pos);
  //   ++pos;
  // }
}

auto beamInterfaceElement2D::getJacXi() -> prec {
  return this->thickness / prec(2);
}

auto beamInterfaceElement2D::getJacEta(ptrCol &pointers, indexType edgeNum)
    -> prec {
  auto &tEdge = pointers.getGeometryData()->getEdgeData(this->edges[edgeNum]);
  auto c1 = tEdge.getCoordinates(-prec(1));
  auto c2 = tEdge.getCoordinates(prec(1));

  auto dir = c2 - c1;
  prec jac = dir.transpose() * this->tangent;
  jac /= prec(2);
  return jac;
}

void beamInterfaceElement2D::getDofsOnVert(ptrCol &pointers,
                                           std::vector<DegreeOfFreedom *> &Dofs,
                                           const indexType &meshID) {

  Geometry::BeamInterface2D *elem;
  elem = this->getInterfaceGeoElem(pointers);
  elem->getDofsOnBeamNode(pointers, Dofs, meshID);
}

void beamInterfaceElement2D::getNodesOnVert(ptrCol &pointers,
                                            std::vector<GenericNodes *> &Nodes,
                                            const indexType &meshID) {
  Nodes.clear();
  auto &tvert = pointers.getGeometryData()->getVertexData(this->vert);
  Nodes = tvert.getNodesOfSet(meshID);
}

auto beamInterfaceElement2D::getNumberOfSolidNodes(ptrCol &pointers,
                                                   const indexType &meshID)
    -> indexType {
  std::vector<GenericNodes *> Nodes;
  this->getNodesOnSolid(pointers, Nodes, meshID);
  auto nn = static_cast<indexType>(Nodes.size());
  return nn;
}

void beamInterfaceElement2D::setBCOnAllNodesSolid(ptrCol &pointers,
                                                  const indexType &meshID,
                                                  const indexType &dof) {

  Geometry::BeamInterface2D *elem;
  elem = this->getInterfaceGeoElem(pointers);
  elem->setBCOnSolid(pointers, meshID, dof);
}

void beamInterfaceElement2D::setBCOnNodeSolid(ptrCol &pointers,
                                              const indexType &meshID,
                                              const indexType &node,
                                              const indexType &dof) {
  std::vector<GenericNodes *> Nodes;
  this->getNodesOnSolid(pointers, Nodes, meshID);
  GenericNodes *tnode;
  tnode = Nodes[node];
  tnode->setBoundaryCondition(dof);
}

void beamInterfaceElement2D::setBCOnVert(ptrCol &pointers,
                                         const indexType &meshID,
                                         const indexType &dof) {

  Geometry::BeamInterface2D *elem;
  elem = this->getInterfaceGeoElem(pointers);
  elem->setBCOnVert(pointers, meshID, dof);
}

void beamInterfaceElement2D::getDofs(ptrCol &pointers,
                                     std::vector<DegreeOfFreedom *> &Dofs,
                                     const indexType &meshIDDisp,
                                     const indexType &meshIDWarp,
                                     const indexType &meshIDRot,
                                     const indexType &order) {
  std::vector<GenericNodes *> tnodes;
  std::vector<GenericNodes *> orderedNodeList;

  indexType counter = 0;
  for (auto &i : this->edges) {
    auto &tedge = pointers.getGeometryData()->getEdgeData(i);
    tedge.getNodes(tnodes, meshIDDisp);
    for (auto &j : tnodes) {
      if (counter == this->dofmapDisp[j->getId()]) {
        ++counter;
        orderedNodeList.push_back(j);
      }
    }
  }
  counter = 0;
  for (auto &i : this->edges) {
    auto &tedge = pointers.getGeometryData()->getEdgeData(i);
    tedge.getNodes(tnodes, meshIDWarp);
    for (auto &j : tnodes) {
      if (counter == this->dofmapWarp[j->getId()]) {
        ++counter;
        orderedNodeList.push_back(j);
      }
    }
  }
  auto &tvert = pointers.getGeometryData()->getVertexData(this->vert);
  tnodes = tvert.getNodesOfSet(meshIDDisp);
  orderedNodeList.push_back(tnodes[0]);
  tnodes = tvert.getNodesOfSet(meshIDRot);
  orderedNodeList.push_back(tnodes[0]);


  for (auto &i : orderedNodeList) {
    auto tdofs = i->getDegreesOfFreedom();
    for (auto &j : tdofs) {
      Dofs.push_back(j);
    }
  }
}

void beamInterfaceElement2D::getDofs2(ptrCol &pointers,
                                      std::vector<DegreeOfFreedom *> &Dofs,
                                      const indexType &meshIDDisp,
                                      const indexType &meshIDRot,
                                      const indexType &order) {
  std::vector<GenericNodes *> tnodes;
  std::vector<GenericNodes *> orderedNodeList;

  indexType counter = 0;
  for (auto &i : this->edges) {
    auto &tedge = pointers.getGeometryData()->getEdgeData(i);
    tedge.getNodes(tnodes, meshIDDisp);
    for (auto &j : tnodes) {
      if (counter == this->dofmapDisp[j->getId()]) {
        ++counter;
        orderedNodeList.push_back(j);
      }
    }
  }

  auto &tvert = pointers.getGeometryData()->getVertexData(this->vert);
  tnodes = tvert.getNodesOfSet(meshIDDisp);
  orderedNodeList.push_back(tnodes[0]);
  tnodes = tvert.getNodesOfSet(meshIDRot);
  orderedNodeList.push_back(tnodes[0]);


  for (auto &i : orderedNodeList) {
    auto tdofs = i->getDegreesOfFreedom();
    for (auto &j : tdofs) {
      Dofs.push_back(j);
    }
  }
}

void beamInterfaceElement2D::setEdges(std::vector<indexType> &edgesIn) {

  for (indexType & i : edgesIn) {
    this->edges.push_back(i);
  }
}

auto beamInterfaceElement2D::getVertexId(ptrCol &pointers, indexType num)
    -> indexType {
  // GeometryData *geom;
  // GenericGeometryElement *geomEle;
  // geom = pointers.getGeometryData();
  // geomEle = geom->getGeometryElement(GeometryTypes::LinearEdge, this->edge);
  // std::vector<indexType> verts;
  // geomEle->getVerts(verts);

  return 0;
}

auto beamInterfaceElement2D::getVertex(ptrCol &pointers, indexType num)
    -> Geometry::VertexData & {

  auto *ele = this->getInterfaceGeoElem(pointers);
  return ele->getVertex(pointers);
}


void beamInterfaceElement2D::computeGeometry(ptrCol &pointers) {
  auto geoData = pointers.getGeometryData();
  auto &tEdge = geoData->getEdgeData(this->edges[0]);

  Types::Vector3<prec> Coordinates1;
  Types::Vector3<prec> Coordinates2;

  Coordinates1 = tEdge.getCoordinates(-1.0);
  Coordinates2 = tEdge.getCoordinates(1.0);

  this->tangent = Coordinates2 - Coordinates1;
  this->tangent /= this->tangent.norm();
  this->normal.setZero();
  this->normal(1) = -this->tangent(0);
  this->normal(0) = this->tangent(1);

  Coordinates2 = geoData->getVertexData(this->vert).getCoordinates();
  this->refCoordinates = Coordinates2 - Coordinates1;
  this->thickness = this->refCoordinates.dot(this->normal);

  this->refCoordinates = Coordinates2 - this->thickness * this->normal;
}

inline void
beamInterfaceElement2D::getH1Dofs(ptrCol &pointers,
                                  std::vector<DegreeOfFreedom *> &Dofs,
                                  indexType meshID, indexType order) {}

void beamInterfaceElement2D::getJacobianInterface(ptrCol &pointers,
                                                  prec &jacobi, const prec &xsi,
                                                  const indexType &edgeNum) {
  auto &tempEdge = pointers.getGeometryData()->getEdgeData(edgeNum);

  Types::Vector3<prec> c1;
  Types::Vector3<prec> c2;
  std::vector<Geometry::GeometryBaseData *> verts;
  c1 = tempEdge.getVertex(0)->getCoordinates();
  c2 = tempEdge.getVertex(1)->getCoordinates();
  jacobi = (c2 - c1).norm() / prec(2);
}

inline void beamInterfaceElement2D::getH1Shapes(
    ptrCol &pointers, indexType order, prec jacobi, Types::VectorX<prec> &shape,
    Types::VectorX<prec> &shapeDerivative, prec xsi) {}

void beamInterfaceElement2D::setShapes(ptrCol &pointers,
                                       const indexType &mesIDDisp,
                                       const indexType &mesIDWarp,
                                       const indexType &mesIDRot,
                                       const indexType &order) {
  for (auto &i : this->edges) {
    auto &tedge = pointers.getGeometryData()->getEdgeData(i);
    tedge.setH1Shapes(mesIDDisp, order, NodeTypes::displacement);
    tedge.setH1Shapes(mesIDWarp, order, NodeTypes::displacement);
  }
  auto tvert = pointers.getGeometryData()->getVertexData(this->vert);
  tvert.setNodeSet(mesIDDisp, 1, NodeTypes::displacement);
  tvert.setNodeSet(mesIDRot, 1, NodeTypes::displacement);

  std::vector<Geometry::GeometryBaseData *> tverts;

  auto &tedge = pointers.getGeometryData()->getEdgeData(this->edges[0]);
  Types::Vector3<prec> c1;
  Types::Vector3<prec> c2;
  Types::Vector3<prec> c3;
  c1 = tedge.getVertex(0)->getCoordinates();
  c2 = tedge.getVertex(1)->getCoordinates();
  c3 = c2 - c1;
  c3 /= c3.norm();
  std::vector<indexType> ind{0, 1};
  this->normal.setZero();
  this->tangent.setZero();
  this->tangent(0) = c3(0);
  this->tangent(1) = c3(1);
  this->normal(0) = c3(1);
  this->normal(1) = -c3(0);

  c1 = pointers.getGeometryData()->getVertexData(this->vert).getCoordinates();
  c1 -= c2;
  this->thickness = c1.dot(this->normal);
  // this->thickness = c1.transpose() * this->normal;
}

void beamInterfaceElement2D::setShapes2(ptrCol &pointers,
                                        const indexType &mesIDDisp,
                                        const indexType &mesIDRot,
                                        const indexType &order) {
  for (auto &i : this->edges) {
    auto &tedge = pointers.getGeometryData()->getEdgeData(i);
    tedge.setH1Shapes(mesIDDisp, order, NodeTypes::displacement);
  }
  auto tvert = pointers.getGeometryData()->getVertexData(this->vert);
  tvert.setNodeSet(mesIDDisp, 1, NodeTypes::displacement);
  tvert.setNodeSet(mesIDRot, 1, NodeTypes::displacement);

  std::vector<Geometry::VertexData *> tverts;

  auto &tedge = pointers.getGeometryData()->getEdgeData(this->edges[0]);
  tedge.getVerts(tverts);
  Types::Vector3<prec> c1;
  Types::Vector3<prec> c2;
  Types::Vector3<prec> c3;
  c1 = tedge.getVertex(0)->getCoordinates(); // tverts[0]->getCoordinates();
  c2 = tedge.getVertex(1)->getCoordinates(); // tverts[1]->getCoordinates();
  c3 = c2 - c1;
  c3 /= c3.norm();
  std::vector<indexType> ind{0, 1};
  this->normal.setZero();
  this->tangent.setZero();
  this->tangent(0) = c3(0);
  this->tangent(1) = c3(1);
  this->normal(0) = c3(1);
  this->normal(1) = -c3(0);

  c1 = pointers.getGeometryData()->getVertexData(this->vert).getCoordinates();
  c1 -= c2;
  this->thickness = c1.dot(this->normal);
  // this->thickness = c1.transpose() * this->normal;
}

void beamInterfaceElement2D::toParaviewAdapter(PointerCollection &ptrCol,
                                               vtkPlotInterface &catalyst,
                                               const ParaviewSwitch &ToDo) {
#ifdef USE_VTK

  switch (ToDo) {
  case ParaviewSwitch::Mesh: {
    Types::Vector3<prec> coors;
    std::vector<double> coorsA;
    coorsA.resize(3);
    std::vector<int> cell(2);


    // Edges to quad
    std::map<indexType, indexType> lIdGIdMap;
    indexType pos = 0;
    for (auto &i : this->edges) {
      auto &tedge = ptrCol.getGeometryData()->getEdgeData(i);
      for (auto vn = 0; vn < tedge.getNumberOfVerts(); ++vn) {
        auto Vert = tedge.getVertex(vn);
        if (lIdGIdMap.find(Vert->getId()) == lIdGIdMap.end()) {
          lIdGIdMap[Vert->getId()] = pos;
          ++pos;
        }
      }
    }
    // indexType nEdgeNodes = lIdGIdMap.size();
    // int cellId = 1;
    // for (auto &i : this->edges) {
    //      tedge = ptrCol.getGeometryData()->getEdge(i);
    //      Geometry::Vertex *VertA;
    //      Geometry::Vertex *VertB;
    //
    //      VertA = tedge->getVertex(ptrCol, 0);
    //      VertB = tedge->getVertex(ptrCol, 1);
    //
    //      indexType idA = lIdGIdMap[VertA->getId()];
    //      indexType idB = lIdGIdMap[VertB->getId()];
    //
    //      indexType idC = idA + nEdgeNodes;
    //      indexType idD = idB + nEdgeNodes;
    //
    //      auto coor = VertA->getCoordinates();
    //      auto tcoor = static_cast<Types::Vector3<prec>>(coor);
    //      catalyst.addPoint(0, matNum, idA, tcoor(0), tcoor(1), tcoor(2));
    //      coor += this->thickness * this->normal;
    //      tcoor = static_cast<Types::Vector3<double>>(coor);
    //      catalyst.addPoint(0, matNum, idC, tcoor(0), tcoor(1), tcoor(2));
    //      coor = VertB->getCoordinates();
    //      tcoor = static_cast<Types::Vector3<double>>(coor);
    //      catalyst.addPoint(0, matNum, idB, tcoor(0), tcoor(1), tcoor(2));
    //      coor += this->thickness * this->normal;
    //      tcoor = static_cast<Types::Vector3<double>>(coor);
    //      catalyst.addPoint(0, matNum, idD, tcoor(0), tcoor(1), tcoor(2));
    //
    //      int celltype = VTK_QUAD;
    //      indexType nn = 4;
    //      indexType elid = (this->id);
    //      std::vector<indexType> cellIds{(idA), (idC), (idD), (idB)};
    //      catalyst.addCell(0, matNum, elid, cellId, cellIds, nn, celltype);
    //      ++cellId;

    // }

    // for (indexType i = 0; i < 2; ++i) {
    //  Vert = this->getVertex(ptrCol, i);
    //  coors = Vert->getCoordinates();
    //  for (auto j = 0; j < 3; ++j) {
    //    coorsA[j] = static_cast<double>(coors(j));
    //  }
    //  int id = static_cast<int>(Vert->getId());
    //  catalyst.addPoint(0, matNum, id, coorsA[0], coorsA[1], coorsA[2]);
    //  cell[i] = static_cast<int>(id);
    //}
    // int celltype = VTK_LINE;
    // int dummy = 1;
    // int nn = 2;
    // int idd = this->id;
    // catalyst.addCell(0, matNum, idd, dummy, &cell[0], &nn, &celltype);
    break;
  }
  case ParaviewSwitch::Solution: {

    int matNum = static_cast<int>(this->getMaterial()->getNumber());

    Types::VectorX<prec> WarpingShapeX;
    Types::VectorX<prec> WarpingShapeY;
    Types::VectorX<prec> WarpingShapeXY;
    Types::VectorX<prec> WarpingShapeYY;

    std::vector<DegreeOfFreedom *> DisplacementDofs;
    std::vector<DegreeOfFreedom *> WarpingDofs;
    std::vector<DegreeOfFreedom *> VertDispDofs;
    std::vector<DegreeOfFreedom *> VertRotDofs;
    

    this->getDofsOnSolid(ptrCol, DisplacementDofs, this->meshIdDisp, 1);
    this->getDofsOnSolid(ptrCol, WarpingDofs, this->meshIdWarp, 1);
    this->getDofsOnVert(ptrCol, VertDispDofs, this->meshIdDisp);
    this->getDofsOnVert(ptrCol, VertRotDofs, this->meshIdRot);


    auto WarpingSolution = ptrCol.getSolutionState()->getSolution(WarpingDofs);
    auto DisplacementSolution =
        ptrCol.getSolutionState()->getSolution(DisplacementDofs);
    auto VertDispSolution =
        ptrCol.getSolutionState()->getSolution(VertDispDofs);
    auto VertRotSolution = ptrCol.getSolutionState()->getSolution(VertRotDofs);

    Types::VectorX<prec> WarpingFunction;
    WarpingFunction.resize(WarpingDofs.size());
    WarpingFunction.setZero();

    // Edges to quad
    std::map<indexType, indexType> lIdGIdMap;
    indexType pos = 0;
    for (auto &i : this->edges) {
      auto &tEdge = ptrCol.getGeometryData()->getEdgeData(i);
      for (auto vn = 0; vn < tEdge.getNumberOfVerts(); ++vn) {
        auto Vert = tEdge.getVertex(vn);
        if (lIdGIdMap.find(Vert->getId()) == lIdGIdMap.end()) {
          lIdGIdMap[Vert->getId()] = pos;
          ++pos;
        }
      }
    }

    for (auto i = 0; i < this->edges.size(); ++i) {
      //auto &tEdge = ptrs.getGeometryData()->getEdge(this->edges[i]);

      this->getWarpingShapes(ptrCol, WarpingShapeX, WarpingShapeY,
                             WarpingShapeXY,
                             WarpingShapeYY, i, -1, 1, this->meshIdWarp);
      WarpingFunction(3 * i) = 0;
      WarpingFunction(3 * i + 1) = 0;
      for (auto j = 0; j < WarpingDofs.size() / 3; ++j) {
        WarpingFunction(3 * i) += WarpingSolution(3 * j) * WarpingShapeX(j);
        WarpingFunction(3 * i + 1) +=
            WarpingSolution(3 * j + 1) * WarpingShapeY(j);
      }

      this->getWarpingShapes(ptrCol, WarpingShapeX, WarpingShapeY,
                             WarpingShapeXY,
                             WarpingShapeYY, i, 1, 1, this->meshIdWarp);
      WarpingFunction(3 * (i + 1)) = 0;
      WarpingFunction(3 * (i + 1) + 1) = 0;
      for (auto j = 0; j < WarpingDofs.size() / 3; ++j) {
        WarpingFunction(3 * (i + 1)) +=
            WarpingSolution(3 * j) * WarpingShapeX(j);
        WarpingFunction(3 * (i + 1) + 1) +=
            WarpingSolution(3 * j + 1) * WarpingShapeY(j);
      }
    }

    std::stringstream ArrName;
    ArrName << "SolutionM" << this->meshIdDisp << "N" << 0;

    indexType nnodes = this->edges.size() + 1;

    std::vector<prec> data;

    Geometry::VertexData *Vert;
    for (auto i = 0; i < this->edges.size(); ++i) {
      auto &tEdge = ptrCol.getGeometryData()->getEdgeData(this->edges[i]);
      prec y = this->getZCoordinate(ptrCol, i, -1);
      Vert = tEdge.getVertex(0);
      indexType parId = lIdGIdMap[Vert->getId()] + nnodes;
      data.clear();
      data.push_back(WarpingFunction(3 * i));
      data.push_back(WarpingFunction(3 * i + 1));
      data.push_back(WarpingFunction(3 * i + 2));
      catalyst.setPointData(0, matNum, parId, data, 3, "Warping Function");
      parId = (lIdGIdMap[Vert->getId()]);
      catalyst.setPointData(0, matNum, parId, data, 3, "Warping Function");

      parId = (lIdGIdMap[Vert->getId()] + nnodes);
      data.clear();
      data.push_back((WarpingFunction(3 * i) + VertDispSolution(0) -
                      y * VertRotSolution(0)));
      data.push_back((WarpingFunction(3 * i + 1) + VertDispSolution(1)));
      data.push_back((WarpingFunction(3 * i + 2)));
      catalyst.setPointData(0, matNum, parId, data, 3, "SolutionM1N0");

      Vert = tEdge.getVertex(1);
      parId = (lIdGIdMap[Vert->getId()] + nnodes);
      data.clear();
      data.push_back((WarpingFunction(3 * (i + 1))));
      data.push_back((WarpingFunction(3 * (i + 1) + 1)));
      data.push_back((WarpingFunction(3 * (i + 1) + 2)));
      catalyst.setPointData(0, matNum, parId, data, 3, "Warping Function");
      parId = (lIdGIdMap[Vert->getId()]);
      catalyst.setPointData(0, matNum, parId, data, 3, "Warping Function");

      parId = static_cast<int>(lIdGIdMap[Vert->getId()] + nnodes);
      data.clear();
      data.push_back((WarpingFunction(3 * (i + 1)) + VertDispSolution(0) -
                      y * VertRotSolution(0)));
      data.push_back((WarpingFunction(3 * (i + 1) + 1) + VertDispSolution(1)));
      data.push_back((WarpingFunction(3 * (i + 1) + 2)));
      catalyst.setPointData(0, matNum, parId, data, 3, "SolutionM1N0");

      Vert = tEdge.getVertex(0);
      parId = (lIdGIdMap[Vert->getId()]);
      data.clear();
      data.push_back((DisplacementSolution(3 * i)));
      data.push_back((DisplacementSolution(3 * i + 1)));
      data.push_back((DisplacementSolution(3 * i + 2)));
      catalyst.setPointData(0, matNum, parId, data, 3, "SolutionM1N0");
      Vert = tEdge.getVertex(1);
      parId = (lIdGIdMap[Vert->getId()]);
      data.clear();
      data.push_back((DisplacementSolution(3 * (i + 1))));
      data.push_back((DisplacementSolution(3 * (i + 1) + 1)));
      data.push_back((DisplacementSolution(3 * (i + 1) + 2)));
      catalyst.setPointData(0, matNum, parId, data, 3, "SolutionM1N0");
    }

    //    int matNum = static_cast<int>(this->getMaterial()->getNumber());
    //    GenericGeometryElement *Vert;
    //    NodeSet *tempSet;
    //    std::vector<DegreeOfFreedom *> Dofs;
    //    Types::VectorX<prec> solution;
    //    std::vector<double> sol(3);
    //    GenericSolutionState *tsol = ptrCol.getSolutionState();
    //    for (indexType i = 0; i < 2; ++i) {
    //      Vert = this->getVertex(ptrCol, i);
    //      std::vector<NodeSet *> sets;
    //      Vert->getSets(ptrCol, sets);
    //      int vertId = static_cast<int>(Vert->getId());
    //      for (auto j = sets.begin(); j != sets.end(); ++j) {
    //        tempSet = *j;
    //        indexType id = tempSet->getMeshId();
    //        indexType nnodes = tempSet->getNumberOfNodes();
    //        for (auto k = 0; k < nnodes; ++k) {
    //          std::stringstream ArrName;
    //          ArrName << "SolutionM" << id << "N" << k;
    //          char *c = new char[ArrName.str().size() + 1];
    //          strcpy(c, ArrName.str().c_str());
    //          Dofs = tempSet->getDegreesOfFreedom(ptrCol);
    //          this->getSolution(ptrCol, Dofs, solution);
    //          for (auto i = 0; i < 3; ++i) {
    //            sol[i] = static_cast<double>(solution(i));
    //          }
    //          int comp = 3;
    //          catalyst.setPointData<double>(0, matNum, vertId, &sol[0], comp,
    //          c); delete[] c;
    //        }
    //      }
    //      if (tsol->getType() >= SolutionTypes::Transient) {
    //        for (auto j = sets.begin(); j != sets.end(); ++j) {
    //          tempSet = *j;
    //          indexType id = tempSet->getMeshId();
    //          indexType nnodes = tempSet->getNumberOfNodes();
    //          for (auto k = 0; k < nnodes; ++k) {
    //            std::stringstream ArrName;
    //            ArrName << "SolutionVelocityM" << id << "N" << k;
    //            char *c = new char[ArrName.str().size() + 1];
    //            strcpy(c, ArrName.str().c_str());
    //            Dofs = tempSet->getDegreesOfFreedom(ptrCol);
    //            this->getVelocity(ptrCol, Dofs, solution);
    //            for (auto i = 0; i < 3; ++i) {
    //              sol[i] = static_cast<double>(solution(i));
    //            }
    //            int comp = 3;
    //            catalyst.setPointData<double>(0, matNum, vertId, &sol[0],
    //            comp, c); delete[] c;
    //          }
    //          for (auto k = 0; k < nnodes; ++k) {
    //            std::stringstream ArrName;
    //            ArrName << "SolutionAccelerationM" << id << "N" << k;
    //            char *c = new char[ArrName.str().size() + 1];
    //            strcpy(c, ArrName.str().c_str());
    //            Dofs = tempSet->getDegreesOfFreedom(ptrCol);
    //            this->getAcceleration(ptrCol, Dofs, solution);
    //            for (auto i = 0; i < 3; ++i) {
    //              sol[i] = static_cast<double>(solution(i));
    //            }
    //            int comp = 3;
    //            catalyst.setPointData<double>(0, matNum, vertId, &sol[0],
    //            comp, c); delete[] c;
    //          }
    //        }
    //      }
    //    }
  } break;
  case ParaviewSwitch::ProjectedValues:
  case ParaviewSwitch::Eigenvector:
  case ParaviewSwitch::Weights: {
  } break;
  }
#endif
}

void beamInterfaceElement2D::setMeshIdDisp(indexType MeshIdDisp) {
  this->meshIdDisp = MeshIdDisp;
}

void beamInterfaceElement2D::setMeshIdWarp(indexType MeshIdWarp) {
  this->meshIdWarp = MeshIdWarp;
}

void beamInterfaceElement2D::setMeshIdRot(indexType MeshIdRot) {
  this->meshIdRot = MeshIdRot;
}

void beamInterfaceElement2D::setSpecial(std::vector<indexType> &specialIn) {
  this->interfaceElementNumber = specialIn[0];
}

void beamInterfaceElement2D::computeWarpingShapesNew(
    PointerCollection &pointers, indexType localOrder) {

  Geometry::BeamInterface2D *elem;

  elem = this->getInterfaceGeoElem(pointers);

  elem->computeWarpingShapes(pointers, localOrder);
}

auto
beamInterfaceElement2D::getInterfaceGeoElem(PointerCollection &pointers) -> Geometry::BeamInterface2D * {

  Geometry::BeamInterface2D *elem;

  elem = static_cast<Geometry::BeamInterface2D *>(
      pointers.getGeometryData()->getSpecial(this->interfaceElementNumber));

  return elem;
}

void beamInterfaceElement2D::getLocalWarpingShapesA1(
    PointerCollection &pointers, Types::VectorX<prec> &shapes,
    Types::VectorX<prec> &shapeDerivative, indexType localEdgeNumber,
    prec eta) {
  Geometry::BeamInterface2D *elem;
  elem = this->getInterfaceGeoElem(pointers);
  elem->getLocalWarpingShapesA1(pointers, shapes,
                                shapeDerivative, localEdgeNumber, eta);
}

void beamInterfaceElement2D::getLocalWarpingShapesA2(
    PointerCollection &pointers, Types::VectorX<prec> &shapes,
    Types::VectorX<prec> &shapeDerivative, indexType localEdgeNumber,
    prec eta) {
  Geometry::BeamInterface2D *elem;
  elem = this->getInterfaceGeoElem(pointers);
  elem->getLocalWarpingShapesA1(pointers, shapes,
                                shapeDerivative, localEdgeNumber, eta);
}

auto
beamInterfaceElement2D::getA1(PointerCollection &pointers) -> Types::Vector3<prec> {
  Geometry::BeamInterface2D *elem;
  elem = this->getInterfaceGeoElem(pointers);

  return elem->getA1();
}

auto
beamInterfaceElement2D::getA2(PointerCollection &pointers) -> Types::Vector3<prec> {
  Geometry::BeamInterface2D *elem;
  elem = this->getInterfaceGeoElem(pointers);

  return elem->getA2();
}

auto beamInterfaceElement2D::getThickness(PointerCollection &pointers) -> prec {
  Geometry::BeamInterface2D *elem;
  elem = this->getInterfaceGeoElem(pointers);

  return elem->getThickness();
}

auto beamInterfaceElement2D::getEdge(PointerCollection &pointers,
                                                 indexType localNumber) -> Geometry::EdgesData & {
  return pointers.getGeometryData()->getEdgeData(this->edges[localNumber]);
}

auto beamInterfaceElement2D::getElementType() -> Elementtypes {
  return Elementtypes::beamInterfaceElement2D;
}

void beamInterfaceElement2D::set_pointers(PointerCollection &pointers) {
}

auto beamInterfaceElement2D::getIntegrationPoints(ptrCol &pointers)
    -> IntegrationPoints {
  auto geoElem = pointers.getGeometryData()->getSpecial(interfaceElementNumber);
  auto GP = geoElem->getIntegrationPoints(this->m_id);
  return GP;
}

} // namespace HierAMuS
