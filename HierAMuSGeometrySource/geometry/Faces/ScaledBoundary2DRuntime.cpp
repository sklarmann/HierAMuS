// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause
#include "geometry/GeometryData.h"

#include "MatrixTypes.h"
#include "geometry/GeometryTypes.h"

#include "plot/vtkplotClass.h"
#include "solver/GenericSolutionState.h"

#include "geometry/Faces/ScaledBoundary2DData.h"
#include "geometry/VertexRuntime.h"
#include "geometry/Edges/EdgesRuntime.h"
#include "geometry/Edges/EdgeH1ShapesInterface.h"
#include <geometry/Faces/ScaledBoundary2DRuntime.h>
#include <sstream>

#include <geometry/Edges/EdgesData.h>
#include <geometry/GeometryData.h>
#include <geometry/VertexData.h>

#include "shapefunctions/IntegrationsPoints/helperClasses/IntegrationPointsManagement.h"

#include "shapefunctions/LobattoShapes.h"

#include <vtkCellType.h>

#include "HelperFunctions.h"

namespace HierAMuS::Geometry {

ScaledBoundary2DRuntime::ScaledBoundary2DRuntime(
    GeometryData &geoData, ScaledBoundary2DData &base_element)
    : FacesRuntime(base_element), m_Scaled_data(base_element){};

ScaledBoundary2DRuntime::~ScaledBoundary2DRuntime() = default;

auto ScaledBoundary2DRuntime::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::ScaledBoundary2DRuntime::type;
}





auto ScaledBoundary2DRuntime::getVertexNumbers() -> std::vector<indexType> {
  return this->verts;
}

void ScaledBoundary2DRuntime::getVertices(
    std::vector<VertexRuntime *> &vertsOut) {
  vertsOut = m_Verts;
}

auto ScaledBoundary2DRuntime::getVertex(indexType local_number)
    -> Geometry::VertexRuntime * {
  return m_Verts[local_number];
}

auto ScaledBoundary2DRuntime::getEdgeNumbers() -> std::vector<indexType> {
  std::vector<indexType> edgesOut;
  edgesOut.insert(edgesOut.end(), edges.begin(), edges.end());
  return edgesOut;
}
void ScaledBoundary2DRuntime::setScalingCenter(prec x, prec y) {
  this->scalingcoor.resize(2);
  this->scalingcoor[0] = x;
  this->scalingcoor[1] = y;
}

void ScaledBoundary2DRuntime::computeScalingCenter() {
  // this->scalingCenter =
  // pointers.getGeometryData()->requestNewGeometryObject(GeometryTypes::Vertex);
  // auto &vertex = pointers.getGeometryData()->getVertex(this->scalingCenter);
  // vertex.setCoordinates(1,2,0);

  // this->scalingcoor.resize(2);
  // this->scalingcoor[0] = 10;
  // this->scalingcoor[1] = 10;

  // ScalingCenter im Kreuz f�r Viereck

  auto tempcoor0 = m_vert_pointers[0]->getCoordinates();
  auto tempcoor1 = m_vert_pointers[1]->getCoordinates();
  auto tempcoor2 = m_vert_pointers[2]->getCoordinates();
  auto tempcoor3 = m_vert_pointers[3]->getCoordinates();

  prec V1x;
  prec V1y;
  prec V2x;
  prec V2y;
  prec V3x;
  prec V3y;
  prec V4x;
  prec V4y;

  V1x = tempcoor0(0);
  V1y = tempcoor0(1);
  V2x = tempcoor1(0);
  V2y = tempcoor1(1);
  V3x = tempcoor2(0);
  V3y = tempcoor2(1);
  V4x = tempcoor3(0);
  V4y = tempcoor3(1);

  prec r;
  prec scx;
  prec scy;
  if (V4y == V2y) {
    r = (V4y - V1y) / (V3y - V1y);
    scx = V1x + r * (V3x - V1x);
    scy = V1y + r * (V3y - V1y);
  } else {

    r = (V1x - V4x + (V3x - V1x) / (V3y - V1y) * (V4y - V1y)) /
        (V2x - V4x - (V3x - V1x) / (V3y - V1y) * (V2y - V4y));
    scx = V1x + r * (V3x - V1x);
    scy = V1y + r * (V3y - V1y);
  }

  this->scalingcoor.resize(2);
  this->scalingcoor[0] = scx;
  this->scalingcoor[1] = scy;

  // Scaling Center im Flaechenschwerpunkt
  /*
  prec A = 0;
  for (auto i = 0; i < (this->edges.size()); ++i) {
      if (i == this->edges.size() - 1) {
          auto tempcoor1 =
  pointers.getGeometryData()->getVertex(this->verts[i]).getCoordinates(); int j
  = i + 1; auto tempcoor2 =
  pointers.getGeometryData()->getVertex(this->verts[0]).getCoordinates(); A +=
  0.5 * (tempcoor1(0) * tempcoor2(1) - tempcoor2(0) * tempcoor1(1));
      }
      else {
          auto tempcoor1 =
  pointers.getGeometryData()->getVertex(this->verts[i]).getCoordinates(); int j
  = i + 1; auto tempcoor2 =
  pointers.getGeometryData()->getVertex(this->verts[j]).getCoordinates(); A +=
  0.5 * (tempcoor1(0) * tempcoor2(1) - tempcoor2(0) * tempcoor1(1));
      }
  }
  //std::cout <<"Fl�che: " <<A<<"\n";
  this->scalingcoor.resize(2);
  for (auto i = 0; i < (this->edges.size()); ++i) {
      if (i == this->edges.size() - 1) {
          auto tempcoor1 =
  pointers.getGeometryData()->getVertex(this->verts[i]).getCoordinates(); int j
  = i + 1; auto tempcoor2 =
  pointers.getGeometryData()->getVertex(this->verts[0]).getCoordinates();
          scalingcoor[0] += ((tempcoor1(0) + tempcoor2(0)) * (tempcoor1(0) *
  tempcoor2(1) - tempcoor2(0) * tempcoor1(1))) / (6 * A); scalingcoor[1] +=
  ((tempcoor1(1) + tempcoor2(1)) * (tempcoor1(0) * tempcoor2(1) - tempcoor2(0) *
  tempcoor1(1))) / (6 * A);
      }
      else {
          auto tempcoor1 =
  pointers.getGeometryData()->getVertex(this->verts[i]).getCoordinates(); int j
  = i + 1; auto tempcoor2 =
  pointers.getGeometryData()->getVertex(this->verts[j]).getCoordinates();
          scalingcoor[0] += ((tempcoor1(0) + tempcoor2(0)) * (tempcoor1(0) *
  tempcoor2(1) - tempcoor2(0) * tempcoor1(1))) / (6 * A); scalingcoor[1] +=
  ((tempcoor1(1) + tempcoor2(1)) * (tempcoor1(0) * tempcoor2(1) - tempcoor2(0) *
  tempcoor1(1))) / (6 * A);
      }
  }
  //std::cout << "Scalingcenter: "<<scalingcoor[0]<<", ";
  //std::cout << scalingcoor[1] << "\n";
  */

  // Scalingcenter aus gemittelten Koordinaten
  /*
  prec x, y, sumx = 0, sumy = 0;
  for (auto i=0; i < this->edges.size();++i) {
      auto tempcoor =
  pointers.getGeometryData()->getVertex(this->verts[i]).getCoordinates(); sumx
  += tempcoor(0); sumy += tempcoor(1);
  }
  x = sumx / this->edges.size();
  y = sumy / this->edges.size();

  this->scalingcoor.resize(2);
  this->scalingcoor[0] = x;
  this->scalingcoor[1] = y;
  */
}

void ScaledBoundary2DRuntime::getEdges(std::vector<EdgesRuntime *> &edgesOut) {
  edgesOut.clear();
  for (auto &edge : this->m_Edges) {
    edgesOut.push_back(edge.get());
  }
}

auto ScaledBoundary2DRuntime::getEdgeNumber(indexType local_number) -> indexType {
  return edges[local_number];
}

auto ScaledBoundary2DRuntime::getEdge(indexType local_number)
    -> Geometry::EdgesRuntime *  {
  return m_Edges[local_number].get();
}

auto ScaledBoundary2DRuntime::hasVertices(indexType v1, indexType v2,
                                          indexType v3) -> bool {
  if (contains(this->verts, v1)) {
    if (contains(this->verts, v2)) {
      if (contains(this->verts, v3)) {
        return true;
      }
    }
  }
  return false;
}

void ScaledBoundary2DRuntime::print(spdlog::logger &Log) {
  m_Scaled_data.print(Log);
}


auto ScaledBoundary2DRuntime::getCoordinates(IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  return {};
}

auto ScaledBoundary2DRuntime::getIntegrationPoints(indexType elementId) -> IntegrationPoints {
  IntegrationPoints temp =
      IntegrationPointsManagement::getIntegrationsPoints(elementId);
  temp.setType(IntegrationType::Scaled2D);
  temp.setNumberOfSections(this->edges.size());
  return temp;
}

auto ScaledBoundary2DRuntime::getJacobian(IntegrationPoint &IntegrationPt) -> Types::Matrix22<prec> {

  Types::Matrix22<prec> jacobi;
  jacobi.setZero();

  auto shapes = this->getH1Shapes(1, IntegrationPt);
  indexType currsection = IntegrationPt.sectionNumber;
  indexType numshapes = shapes.shapes.rows();
  jacobi(0, 0) += shapes.shapeDeriv(0, numshapes - 1) * this->scalingcoor[0];
  jacobi(0, 1) += shapes.shapeDeriv(1, numshapes - 1) * this->scalingcoor[0];
  jacobi(1, 0) += shapes.shapeDeriv(0, numshapes - 1) * this->scalingcoor[1];
  jacobi(1, 1) += shapes.shapeDeriv(1, numshapes - 1) * this->scalingcoor[1];

  if ((currsection + 1) == edges.size()) {
    // int a = currsection + 1;
    auto coord = m_Verts[currsection]->getCoordinates();
    jacobi(0, 0) += shapes.shapeDeriv(0, currsection) * coord(0);
    jacobi(0, 1) += shapes.shapeDeriv(1, currsection) * coord(0);
    jacobi(1, 0) += shapes.shapeDeriv(0, currsection) * coord(1);
    jacobi(1, 1) += shapes.shapeDeriv(1, currsection) * coord(1);

    coord = m_Verts[0]->getCoordinates();
    jacobi(0, 0) += shapes.shapeDeriv(0, 0) * coord(0);
    jacobi(0, 1) += shapes.shapeDeriv(1, 0) * coord(0);
    jacobi(1, 0) += shapes.shapeDeriv(0, 0) * coord(1);
    jacobi(1, 1) += shapes.shapeDeriv(1, 0) * coord(1);
  } else {
    for (auto i = currsection; i < currsection + 2; ++i) {

      auto coord = m_Verts[i]->getCoordinates();
      jacobi(0, 0) += shapes.shapeDeriv(0, i) * coord(0);
      jacobi(0, 1) += shapes.shapeDeriv(1, i) * coord(0);
      jacobi(1, 0) += shapes.shapeDeriv(0, i) * coord(1);
      jacobi(1, 1) += shapes.shapeDeriv(1, i) * coord(1);
    }
  }
  return jacobi;
}

// H1Shapes

void ScaledBoundary2DRuntime::setH1Shapes(indexType meshId, indexType order,
                                          NodeTypes type) {
  for (auto edge : this->edges) {
    auto &edgeTemp = *m_Edges[edge]->getH1Edge();
    edgeTemp.setH1Shapes(meshId, order, type);
  }
  this->setH1ShapesInternal(meshId, order, type);
}

void ScaledBoundary2DRuntime::setH1ShapesInternal(indexType meshId,
                                                  indexType order,
                                                  NodeTypes type) {
  // ScalingCenter
  if (order == 1) {
    this->setNodeSet(meshId, 1, type);
  }

  if (order > 1) {
    indexType numNodes = 1 + (order - 1) * (this->verts.size() +
                                            (order - 1) * this->edges.size());

    this->setNodeSet(meshId, numNodes, type);
  }
}

void ScaledBoundary2DRuntime::getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs,
                                        indexType meshID, indexType order) {

  for (auto i = 0; i < this->verts.size(); ++i) {
    auto &tempVert = *m_Verts[i];
    auto nodeList = tempVert.getNodeSetNodeListMeshId(meshID);
    auto tdofs(nodeList.getDegreesOfFreedom());
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
  this->getH1DofsInternal(Dofs, meshID, order);
}

auto ScaledBoundary2DRuntime::getH1Dofs(indexType meshID, indexType order)
    -> std::vector<DegreeOfFreedom *>  {
  std::vector<DegreeOfFreedom *> Dofs;
  this->getH1Dofs(Dofs, meshID, order);
  return Dofs;
}

void ScaledBoundary2DRuntime::getH1DofsInternal(
    std::vector<DegreeOfFreedom *> &Dofs,
    indexType meshID, indexType order) {
  // ScalingCenter
  if (order == 1) {
    auto nodeList = this->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }
  if (order > 1) {
    for (auto edge : this->edges) {
      auto &tempEdge = *m_Edges[edge]->getH1Edge();
      tempEdge.getH1DofsInternal(Dofs, meshID, order);
    }
    // innere Nodes
    auto nodeList = this->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }
}

auto ScaledBoundary2DRuntime::getH1NodesList(indexType meshID, indexType order)
    -> MeshIdNodeList {
  return MeshIdNodeList(meshID);
}

auto ScaledBoundary2DRuntime::getH1ShapesInternal(
    indexType order, IntegrationPoint &ip,
    faceorientation orientation) -> H1Shapes {

  if (order > 1) {
    indexType numnodes = (order + 1) * (order - 1);
    auto shapes = H1Shapes(numnodes, 2);

    prec ni;
    prec dni;
    prec sj;
    prec dsj;
    indexType counter;
    indexType rowcount;
    indexType gradxsi;
    indexType gradeta;
    rowcount = 0;

    prec xsi = ip.xi;
    prec eta = ip.eta;

    for (auto i = 0; i < order - 1; ++i) {
      gradxsi = i + 2;
      LobattoShapes::getShape(ni, dni, xsi, gradxsi);

      // linear
      LobattoShapes::getShape(sj, dsj, eta, 0);
      shapes.shapes(i + rowcount) = ni * sj;
      shapes.shapeDeriv(0, i + rowcount) = sj * dni;
      shapes.shapeDeriv(1, i + rowcount) = ni * dsj;
      LobattoShapes::getShape(sj, dsj, eta, 1);
      shapes.shapes(i + rowcount + order) = ni * sj;
      shapes.shapeDeriv(0, i + rowcount + order) = dni * sj;
      shapes.shapeDeriv(1, i + rowcount + order) = ni * dsj;

      // higher funkt
      for (auto j = 0; j < order - 1; ++j) {
        gradeta = j + 2;
        LobattoShapes::getShape(sj, dsj, eta, gradeta);
        counter = 1 + i + j + rowcount;
        shapes.shapes(counter) = ni * sj;
        shapes.shapeDeriv(0, counter) = dni * sj;
        shapes.shapeDeriv(1, counter) = ni * dsj;
      }

      rowcount += order;
    }
    return shapes;
  } else {
    return H1Shapes(0, 2);
  }
}

auto ScaledBoundary2DRuntime::getH1Nodes(indexType meshID, indexType order)
    -> std::vector<GenericNodes *>  {
  return std::vector<GenericNodes *>();
}

auto ScaledBoundary2DRuntime::getH1NodesInternal(indexType meshID,
                                                 indexType order)
    -> std::vector<GenericNodes *> {
  return std::vector<GenericNodes *>();
}

auto ScaledBoundary2DRuntime::getH1Shapes(indexType order,
                                          IntegrationPoint &IntegrationPt) -> H1Shapes {


  indexType nRand = order * this->edges.size();
  indexType nInnen = nRand * (order - 1);
  indexType numshapes = nRand + nInnen + 1;

  H1Shapes shapes(numshapes, 2);

  prec n1;
  prec dn1;
  prec n2;
  prec dn2;
  prec s1;
  prec ds1;
  prec s2;
  prec ds2;
  prec xsi = IntegrationPt.xi;
  prec eta = IntegrationPt.eta;

  if ((IntegrationPt.sectionNumber + 1) != this->edges.size()) {
    // linear

    LobattoShapes::getShape(n1, dn1, xsi, 0);
    LobattoShapes::getShape(n2, dn2, xsi, 1);
    LobattoShapes::getShape(s1, ds1, eta, 0);
    LobattoShapes::getShape(s2, ds2, eta, 1);

    shapes.shapes(numshapes - 1) = n1;
    shapes.shapes(IntegrationPt.sectionNumber) = n2 * s1;
    shapes.shapes(IntegrationPt.sectionNumber + 1) = n2 * s2;

    shapes.shapeDeriv(0, numshapes - 1) = dn1;
    shapes.shapeDeriv(0, IntegrationPt.sectionNumber) = dn2 * s1;
    shapes.shapeDeriv(0, IntegrationPt.sectionNumber + 1) = dn2 * s2;

    shapes.shapeDeriv(1, numshapes - 1) = 0;
    shapes.shapeDeriv(1, IntegrationPt.sectionNumber) = n2 * ds1;
    shapes.shapeDeriv(1, IntegrationPt.sectionNumber + 1) = n2 * ds2;

    // addition weiterer funktionen
    if (order > 1) {
      Types::VectorX<prec> tempshape;
      Types::VectorX<prec> tempEdgeDerivative;
      Types::Matrix2X<prec> tempshapeDerivative;
      // Edge
      indexType counter =
          this->edges.size() + IntegrationPt.sectionNumber * (order - 1);
      {
        IntegrationPoint eint;
        eint.xi = IntegrationPt.eta;
        auto eshapes = 
        m_Edges[IntegrationPt.sectionNumber]
                           ->getH1Edge()
                           ->getH1ShapesInternal(order, eint);

        for (auto i = 0; i < tempshape.rows(); ++i) {
          shapes.shapes(counter) = eshapes.shapes(i) * n2;
          shapes.shapeDeriv(0, counter) = eshapes.shapes(i) * dn2;
          shapes.shapeDeriv(1, counter) = eshapes.shapeDeriv(i) * n2;
          ++counter;
        }
      }
      // Internal
      auto ishapes = this->getH1ShapesInternal(order, IntegrationPt);
      counter = nRand + IntegrationPt.sectionNumber * order;

      for (auto i = 0; i < tempshape.rows(); ++i) {

        shapes.shapes(counter) = tempshape(i);
        shapes.shapeDeriv(0, counter) = tempshapeDerivative(0, i);
        shapes.shapeDeriv(1, counter) = tempshapeDerivative(1, i);
        if ((i + 1) % (order + 1) == 0) {
          counter += nRand - order;
        } else {
          ++counter;
        }
      }
    }
  } else {
    // linear

    LobattoShapes::getShape(n1, dn1, xsi, 0);
    LobattoShapes::getShape(n2, dn2, xsi, 1);
    LobattoShapes::getShape(s1, ds1, eta, 0);
    LobattoShapes::getShape(s2, ds2, eta, 1);

    shapes.shapes(numshapes - 1) = n1;
    shapes.shapes(IntegrationPt.sectionNumber) = n2 * s1;
    shapes.shapes(0) = n2 * s2;

    shapes.shapeDeriv(0, numshapes - 1) = dn1;
    shapes.shapeDeriv(0, IntegrationPt.sectionNumber) = dn2 * s1;
    shapes.shapeDeriv(0, 0) = dn2 * s2;

    shapes.shapeDeriv(1, numshapes - 1) = 0;
    shapes.shapeDeriv(1, IntegrationPt.sectionNumber) = n2 * ds1;
    shapes.shapeDeriv(1, 0) = n2 * ds2;

    // addition of further functions
    if (order > 1) {
      // Edge
      indexType counter =
          this->edges.size() + IntegrationPt.sectionNumber * (order - 1);
      {
        IntegrationPoint eint;
        eint.xi = IntegrationPt.eta;
        auto eshape = m_Edges[IntegrationPt.sectionNumber]
                          ->getH1Edge()
                          ->getH1ShapesInternal(order, eint);

        for (auto i = 0; i < eshape.shapes.rows(); ++i) {
          shapes.shapes(counter) = eshape.shapes(i) * n2;
          shapes.shapeDeriv(0, counter) = eshape.shapes(i) * dn2;
          shapes.shapeDeriv(1, counter) = eshape.shapeDeriv(i) * n2;
          ++counter;
        }
      }
      // Internal
      auto ishape = this->getH1ShapesInternal(order, IntegrationPt);
      counter = nRand + IntegrationPt.sectionNumber * order;
      indexType rowcounter = 0;
      for (auto i = 0; i < ishape.shapes.rows(); ++i) {
        if ((i + 1) % (order + 1) == 0) {
          shapes.shapes(nRand + rowcounter) = ishape.shapes(i);
          shapes.shapeDeriv(0, nRand + rowcounter) = ishape.shapeDeriv(0, i);
          shapes.shapeDeriv(1, nRand + rowcounter) = ishape.shapeDeriv(1, i);
          rowcounter += nRand;
          counter += nRand - order;
        } else {
          shapes.shapes(counter) = ishape.shapes(i);
          shapes.shapeDeriv(0, counter) = ishape.shapeDeriv(0, i);
          shapes.shapeDeriv(1, counter) = ishape.shapeDeriv(1, i);
          ++counter;
        }
      }
    }
  }

  return shapes;
}



// L2 Shapes

void ScaledBoundary2DRuntime::getL2Dofs(std::vector<DegreeOfFreedom *> &Dofs,
                                        indexType meshID, indexType order) {
  auto nodeList = this->getNodeSetNodeListMeshId(meshID);
  auto tdofs = nodeList.getDegreesOfFreedom();
  Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
}

void ScaledBoundary2DRuntime::setL2Shapes(indexType meshId, indexType order,
                                          NodeTypes type) {
  if (order == 0) {
    this->setNodeSet(meshId, 1, type);
  }
  // nur Knoten an Kante und im Inneren,Knoten an Vertexe nicht beachtet
  if (order >= 1) {
    // indexType numNodes = 1 + (order - 1) * this->edges.size() + (order - 1) *
    // (this->verts.size() + (order - 1) * this->edges.size());
    indexType numNodes =
        (this->verts.size() + (order - 1) * this->edges.size()) * order + 1;
    this->setNodeSet(meshId, numNodes, type);
  }
}

auto ScaledBoundary2DRuntime::getL2Shapes(indexType order,
                                          IntegrationPoint &IntegrationPt)
    -> L2Shapes {
  L2Shapes shapes;

  indexType nRand = order * this->edges.size();
  indexType nInnen = nRand * (order - 1);
  indexType numshapes = nRand + nInnen + 1;

  shapes.shapes.resize(numshapes);
  shapes.shapes.setZero();

  prec n1;
  prec dn1;
  prec n2;
  prec dn2;
  prec s1;
  prec ds1;
  prec s2;
  prec ds2;
  prec xsi = IntegrationPt.xi;
  prec eta = IntegrationPt.eta;

  if ((IntegrationPt.sectionNumber + 1) != this->edges.size()) {
    // linear

    LobattoShapes::getShape(n1, dn1, xsi, 0);
    LobattoShapes::getShape(n2, dn2, xsi, 1);
    LobattoShapes::getShape(s1, ds1, eta, 0);
    LobattoShapes::getShape(s2, ds2, eta, 1);

    shapes.shapes(numshapes - 1) = n1;
    shapes.shapes(IntegrationPt.sectionNumber) = n2 * s1;
    shapes.shapes(IntegrationPt.sectionNumber + 1) = n2 * s2;

    // addition weiterer funktionen
    if (order > 1) {
      Types::VectorX<prec> tempshape;
      Types::VectorX<prec> tempEdgeDerivative;
      Types::Matrix2X<prec> tempshapeDerivative;
      prec si;
      prec dsi;
      // Edge
      indexType counter =
          this->edges.size() + IntegrationPt.sectionNumber * (order - 1);
      {
        tempshape.resize(order - 1);
        for (auto i = 0; i < (order - 1); ++i) {
          indexType gradeta = i + 2;
          LobattoShapes::getShape(si, dsi, eta, gradeta);
          tempshape(i) = si;
        }

        for (auto i = 0; i < tempshape.rows(); ++i) {
          shapes.shapes(counter) = tempshape(i) * n2;
          ++counter;
        }
      }
      // Internal
      this->getL2ShapesInternal(order, tempshape, tempshapeDerivative,
                                xsi, eta);
      counter = nRand + IntegrationPt.sectionNumber * order;

      for (auto i = 0; i < tempshape.rows(); ++i) {

        shapes.shapes(counter) = tempshape(i);
        if ((i + 1) % (order + 1) == 0) {
          counter += nRand - order;
        } else {
          ++counter;
        }
      }
    }
  } else {
    // linear

    LobattoShapes::getShape(n1, dn1, xsi, 0);
    LobattoShapes::getShape(n2, dn2, xsi, 1);
    LobattoShapes::getShape(s1, ds1, eta, 0);
    LobattoShapes::getShape(s2, ds2, eta, 1);

    shapes.shapes(numshapes - 1) = n1;
    shapes.shapes(IntegrationPt.sectionNumber) = n2 * s1;
    shapes.shapes(0) = n2 * s2;

    // addition weiterer funktionen
    if (order > 1) {
      Types::VectorX<prec> tempshape;
      Types::VectorX<prec> tempEdgeDerivative;
      Types::Matrix2X<prec> tempshapeDerivative;
      prec si;
      prec dsi;
      // Edge
      indexType counter =
          this->edges.size() + IntegrationPt.sectionNumber * (order - 1);
      {
        tempshape.resize(order - 1);
        for (auto i = 0; i < (order - 1); ++i) {
          indexType gradeta = i + 2;
          LobattoShapes::getShape(si, dsi, eta, gradeta);
          tempshape(i) = si;
        }

        for (auto i = 0; i < tempshape.rows(); ++i) {
          shapes.shapes(counter) = tempshape(i) * n2;
          ++counter;
        }
      }
      // Internal
      this->getL2ShapesInternal(order, tempshape, tempshapeDerivative,
                                xsi, eta);
      counter = nRand + IntegrationPt.sectionNumber * order;
      indexType rowcounter = 0;
      for (auto i = 0; i < tempshape.rows(); ++i) {
        if ((i + 1) % (order + 1) == 0) {
          shapes.shapes(nRand + rowcounter) = tempshape(i);
          rowcounter += nRand;
          counter += nRand - order;
        } else {
          shapes.shapes(counter) = tempshape(i);
          ++counter;
        }
      }
    }
  }

  return shapes;
}

void ScaledBoundary2DRuntime::getL2ShapesInternal(
    indexType order, Types::VectorX<prec> &shape,
    Types::Matrix2X<prec> &shapeDerivative, prec xsi, prec eta) {

  if (order > 1) {
    indexType numnodes = (order + 1) * (order - 1);
    shape.resize(numnodes);
    shapeDerivative.resize(2, numnodes);

    prec ni;
    prec dni;
    prec sj;
    prec dsj;
    indexType counter;
    indexType rowcount;
    indexType gradxsi;
    indexType gradeta;
    rowcount = 0;

    for (auto i = 0; i < order - 1; ++i) {
      gradxsi = i + 2;
      LobattoShapes::getShape(ni, dni, xsi, gradxsi);

      // linear
      LobattoShapes::getShape(sj, dsj, eta, 0);
      shape(i + rowcount) = ni * sj;
      LobattoShapes::getShape(sj, dsj, eta, 1);
      shape(i + rowcount + order) = ni * sj;

      // higher funkt
      for (auto j = 0; j < order - 1; ++j) {
        gradeta = j + 2;
        LobattoShapes::getShape(sj, dsj, eta, gradeta);
        counter = 1 + i + j + rowcount;
        shape(counter) = ni * sj;
      }
      rowcount += order;
    }
  } else {
    shape.resize(0);
    shapeDerivative.resize(2, 0);
  }
}

void ScaledBoundary2DRuntime::setAllNodeBoundaryConditionMeshId(
    indexType meshId, indexType dof) {
  m_Scaled_data.setAllNodeBoundaryConditionMeshId(meshId, dof);
}

void ScaledBoundary2DRuntime::geometryToParaview(
    vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh) {
  indexType numPoints = this->edges.size();
  std::vector<indexType> points(numPoints);
  points.clear();
  for (auto i : this->verts) {
    auto &vert = *m_vert_pointers[i];
    vert.geometryToParaview(paraviewAdapter, mainMesh, subMesh);
    points.push_back(i);
  }

  paraviewAdapter.addCell(mainMesh, subMesh, this->m_Scaled_data.getId(), 1,
                          points, numPoints, VTK_POLYGON);
}

void ScaledBoundary2DRuntime::computeWeightsParaview(
    vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh) {

  auto GP = this->getIntegrationPoints(-1);
  GP.setOrder(2);

  for (auto i : GP) {
    auto jaco = this->getJacobian(i);
    auto shapes = this->getH1Shapes(1, i);
    prec dA = jaco.determinant() * i.weight;

    for (auto i = 0; i < 4; ++i) {
      auto &vert = *m_vert_pointers[i];
      std::vector<prec> val;
      val.push_back(shapes.shapes(i) * dA);
      paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, val, vert.getId(),
                                           1, paraviewNames::weightName());
    }
  }
}

void ScaledBoundary2DRuntime::H1SolutionToParaview(
    vtkPlotInterface &paraviewAdapter, indexType mainMesh, indexType subMesh,
    indexType order, Types::VectorX<prec> &solution, std::string &name) {
  for (auto i = 0; i < this->edges.size(); ++i) {
    auto &V = *m_vert_pointers[i];
    std::vector<prec> sol(3);
    for (auto j = 0; j < 3; ++j) {
      sol[j] = solution(3 * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, V.getId(), sol, 3, name);
  }
}
void ScaledBoundary2DRuntime::H1DataToParaview(
    vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, Types::VectorX<prec> &Data,
    indexType numberComponents, indexType order, std::string &name) {
  for (auto i = 0; i < this->edges.size(); ++i) {
    auto &V = *m_vert_pointers[i];
    std::vector<prec> sol(numberComponents);
    for (auto j = 0; j < numberComponents; ++j) {
      sol[j] = Data(numberComponents * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, V.getId(), sol,
                                 numberComponents, name);
  }
}

void ScaledBoundary2DRuntime::projectDataToParaviewVertices(
    vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {

  auto shapes = this->getH1Shapes(1, IntegrationPt);
  std::vector<prec> vals(numberComponents);

  auto jaco = this->getJacobian(IntegrationPt);
  auto dA = jaco.determinant() * IntegrationPt.weight;
  for (auto i = 0; i < this->edges.size(); ++i) {
    auto &V = *m_vert_pointers[i];
    for (auto j = 0; j < numberComponents; ++j) {
      vals[j] = data(j) * shapes.shapes(i) * dA;
    }
    paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, vals, V.getId(),
                                         numberComponents, name);
  }
}

void ScaledBoundary2DRuntime::flip() {
  throw std::runtime_error(
      "ERROR in ScaledBoundary2D::flip: Method not implemented!");
}

void ScaledBoundary2DRuntime::rotate(indexType n) {
  throw std::runtime_error(
      "ERROR in ScaledBoundary2D::rotate: Method not implemented!");
}

auto ScaledBoundary2DRuntime::computeMeanCoordinate()
    -> Types::Vector3<prec> {
  throw std::runtime_error("ERROR in ScaledBoundary2D::computeMeanCoordinate: "
                           "Method not implemented!");
  return Types::Vector3<prec>();
}

auto ScaledBoundary2DRuntime::getFaceNormal()
    -> Types::Vector3<prec> {
  throw std::runtime_error("ERROR in ScaledBoundary2D::getFaceNormal: "
                           "Method not implemented!");
  return Types::Vector3<prec>();
}

void ScaledBoundary2DRuntime::set_geometry_pointers(GeometryData &geoData) {
  m_edge_pointers.clear();
  m_vert_pointers.clear();
  for (auto i : verts) {
    m_vert_pointers.push_back(&geoData.getVertexData(i));
  }
  for (auto i : edges) {
    m_edge_pointers.push_back(&geoData.getEdgeData(i));
  }
}

auto ScaledBoundary2DRuntime::getVertexNumber(indexType localNumber)
    -> indexType {
  return verts[localNumber];
}

const GeometryTypes ScaledBoundary2DRuntime::type =
    GeometryTypes::ScaledBoundary2D;
} // namespace HierAMuS::Geometry
