// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "geometry/GeometryData.h"
#include "MatrixTypes.h"
#include "geometry/GeometryTypes.h"

#include "plot/vtkplotClass.h"

#include "geometry/Faces/ScaledBoundary2DRuntime.h"
#include <geometry/Faces/ScaledBoundary2DData.h>
#include <sstream>

#include <geometry/Edges/EdgesData.h>
#include "geometry/GeometryData.h"
#include <geometry/VertexData.h>

#include "shapefunctions/LobattoShapes.h"

#include <vtkCellType.h>

#include "HelperFunctions.h"
#include "shapefunctions/IntegrationsPoints/helperClasses/IntegrationPointsManagement.h"

namespace HierAMuS::Geometry {

ScaledBoundary2DData::ScaledBoundary2DData() = default;

ScaledBoundary2DData::~ScaledBoundary2DData() = default;

auto ScaledBoundary2DData::getRuntimeObject(GeometryData &geoData)
    -> std::shared_ptr<FacesRuntime> {
  return std::make_shared<ScaledBoundary2DRuntime>(geoData, *this);
}

auto ScaledBoundary2DData::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::ScaledBoundary2DData::type;
}

void ScaledBoundary2DData::setVerts(GeometryData &geoData,
                                    std::vector<indexType> &vertsIn) {
  indexType numVerts = vertsIn.size();
  if (numVerts >= 3) {
    this->m_verts.resize(numVerts);
    for (auto i = 0; i < vertsIn.size(); ++i) {
      this->m_verts[i] = vertsIn[i];
    }
  } else {
    std::stringstream temp;
    temp
        << "ScaledBoundary2D geometry element needs atleast 3 vertices. Cannot "
        << "add the given amount of " << vertsIn.size() << " vertices!";
    std::runtime_error(temp.str());
  }
}

void ScaledBoundary2DData::setEdges(const std::vector<indexType> &edgesIn) {
  indexType numEdges = edgesIn.size();
  if (edgesIn.size() >= 3) {
    this->m_edges.resize(numEdges);
    for (auto i = 0; i < edgesIn.size(); ++i) {
      this->m_edges[i] = edgesIn[i];
    }
  } else {
    std::stringstream temp;
    temp << "ScaledBoundary2D geometry element needs at least 3 vertices. "
            "Cannot "
         << "add the given amount of " << edgesIn.size() << " edges!";
    std::runtime_error(temp.str());
  }
}

auto ScaledBoundary2DData::getVertexNumbers() -> std::vector<indexType> {
  return this->m_verts;
}

void ScaledBoundary2DData::getVerts(std::vector<VertexData *> &vertsOut) {
  vertsOut =
      std::vector<VertexData *>(m_vert_pointers.begin(), m_vert_pointers.end());
}

auto ScaledBoundary2DData::getEdgeNumbers() -> std::vector<indexType> {
  std::vector<indexType> edgesOut;
  edgesOut.insert(edgesOut.end(), m_edges.begin(), m_edges.end());
  return edgesOut;
}
void ScaledBoundary2DData::setScalingCenter(prec x, prec y) {
  this->scalingcoor.resize(2);
  this->scalingcoor[0] = x;
  this->scalingcoor[1] = y;
}

void ScaledBoundary2DData::computeScalingCenter() {
  // this->scalingCenter =
  // pointers.getGeometryData()->requestNewGeometryObject(GeometryTypes::Vertex);
  // auto &vertex = pointers.getGeometryData()->getVertex(this->scalingCenter);
  // vertex.setCoordinates(1,2,0);

  // this->scalingcoor.resize(2);
  // this->scalingcoor[0] = 10;
  // this->scalingcoor[1] = 10;

  // ScalingCenter im Kreuz f�r Viereck

  auto tempcoor0 = this->m_vert_pointers[0]->getCoordinates();
  auto tempcoor1 = this->m_vert_pointers[1]->getCoordinates();
  auto tempcoor2 = this->m_vert_pointers[2]->getCoordinates();
  auto tempcoor3 = this->m_vert_pointers[3]->getCoordinates();

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

void ScaledBoundary2DData::getEdges(std::vector<EdgesData *> &edgesOut) {
  edgesOut = std::vector<EdgesData *>(m_edges_pointers.begin(),
                                      m_edges_pointers.end());
}

auto ScaledBoundary2DData::getEdge(indexType local_number) -> EdgesData * {
  return m_edges_pointers[local_number];
}

auto ScaledBoundary2DData::hasVertices(indexType v1, indexType v2, indexType v3)
    -> bool {
  if (contains(this->m_verts, v1)) {
    if (contains(this->m_verts, v2)) {
      if (contains(this->m_verts, v3)) {
        return true;
      }
    }
  }
  return false;
}

void ScaledBoundary2DData::print(spdlog::logger &Logger) {
  Logger.debug("Scaled Boundary 2D Face id:   {:>10}", this->id);
  Logger.debug("Vertices:  {}", fmt::join(this->m_verts, " "));
  Logger.debug("Edges:     {}", fmt::join(this->m_edges, " "));
  this->printEqInfo(Logger);
}

auto ScaledBoundary2DData::getCoordinates(prec xi, prec eta)
    -> Types::Vector3<prec> {
  return {};
}

auto ScaledBoundary2DData::getCoordinates(IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  return {};
}

auto ScaledBoundary2DData::getIntegrationPoints(indexType elementId)
    -> IntegrationPoints {
  IntegrationPoints temp =
      IntegrationPointsManagement::getIntegrationsPoints(elementId);
  temp.setType(IntegrationType::Scaled2D);
  temp.setNumberOfSections(this->m_edges.size());
  return temp;
}

auto ScaledBoundary2DData::getJacobian(IntegrationPoint &IntegrationPt)
    -> Types::Matrix22<prec> {

  Types::Matrix22<prec> jacobi;
  jacobi.setZero();

  auto shapes = this->getH1Shapes(1, IntegrationPt);
  indexType currsection = IntegrationPt.sectionNumber;
  indexType numshapes = shapes.shapes.rows();
  jacobi(0, 0) += shapes.shapeDeriv(0, numshapes - 1) * this->scalingcoor[0];
  jacobi(0, 1) += shapes.shapeDeriv(1, numshapes - 1) * this->scalingcoor[0];
  jacobi(1, 0) += shapes.shapeDeriv(0, numshapes - 1) * this->scalingcoor[1];
  jacobi(1, 1) += shapes.shapeDeriv(1, numshapes - 1) * this->scalingcoor[1];

  if ((currsection + 1) == m_edges.size()) {
    // int a = currsection + 1;
    auto coord = m_vert_pointers[currsection]->getCoordinates();
    jacobi(0, 0) += shapes.shapeDeriv(0, currsection) * coord(0);
    jacobi(0, 1) += shapes.shapeDeriv(1, currsection) * coord(0);
    jacobi(1, 0) += shapes.shapeDeriv(0, currsection) * coord(1);
    jacobi(1, 1) += shapes.shapeDeriv(1, currsection) * coord(1);

    coord = m_vert_pointers[0]->getCoordinates();
    jacobi(0, 0) += shapes.shapeDeriv(0, 0) * coord(0);
    jacobi(0, 1) += shapes.shapeDeriv(1, 0) * coord(0);
    jacobi(1, 0) += shapes.shapeDeriv(0, 0) * coord(1);
    jacobi(1, 1) += shapes.shapeDeriv(1, 0) * coord(1);
  } else {
    for (auto i = currsection; i < currsection + 2; ++i) {

      auto coord = m_vert_pointers[i]->getCoordinates();

      jacobi(0, 0) += shapes.shapeDeriv(0, i) * coord(0);
      jacobi(0, 1) += shapes.shapeDeriv(1, i) * coord(0);
      jacobi(1, 0) += shapes.shapeDeriv(0, i) * coord(1);
      jacobi(1, 1) += shapes.shapeDeriv(1, i) * coord(1);
    }
  }
  return jacobi;
}

// H1Shapes

void ScaledBoundary2DData::setH1Shapes(indexType meshId, indexType order,
                                       NodeTypes type) {
  for (auto edge : this->m_edges) {
    auto &edgeTemp = *m_edges_pointers[edge];
    edgeTemp.setH1Shapes(meshId, order, type);
  }
  this->setH1ShapesInternal(meshId, order, type);
}

void ScaledBoundary2DData::setH1ShapesInternal(indexType meshId,
                                               indexType order,
                                               NodeTypes type) {
  // ScalingCenter
  if (order == 1) {
    this->setNodeSet(meshId, 1, type);
  }

  if (order > 1) {
    indexType numNodes = 1 + (order - 1) * (this->m_verts.size() +
                                            (order - 1) * this->m_edges.size());

    this->setNodeSet(meshId, numNodes, type);
  }
}

void ScaledBoundary2DData::getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs,
                                     indexType meshID, indexType order) {

  for (auto i = 0; i < this->m_verts.size(); ++i) {
    auto nodeList = m_vert_pointers[i]->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }
  this->getH1DofsInternal(Dofs, meshID, order);
}

void ScaledBoundary2DData::getH1DofsInternal(
    std::vector<DegreeOfFreedom *> &Dofs, indexType meshID, indexType order) {
  // ScalingCenter
  if (order == 1) {
    auto nodeList = this->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }
  if (order > 1) {
    for (auto edge : this->m_edges_pointers) {
      edge->getH1DofsInternal(Dofs, meshID, order);
    }
    // innere Nodes
    auto nodeList = this->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }
}

auto ScaledBoundary2DData::getH1NodesList(indexType meshID, indexType order)
    -> MeshIdNodeList {
  return MeshIdNodeList(meshID);
}

void ScaledBoundary2DData::getH1Shapes(indexType order,
                                       Types::VectorX<prec> &shape,
                                       Types::Matrix2X<prec> &shapeDerivative,
                                       prec xsi, prec eta) {
  // fehlerbehaftet

  indexType nRand = order * this->m_edges.size();
  indexType nInnen = nRand * (order - 1);
  indexType numshapes = nRand + nInnen + 1;

  shape.resize(numshapes);
  shape.setZero();
  shapeDerivative.resize(2, numshapes);
  shapeDerivative.setZero();

  prec n1;
  prec dn1;
  prec n2;
  prec dn2;
  prec s1;
  prec ds1;
  prec s2;
  prec ds2;

  // linear
  LobattoShapes::getShape(n1, dn1, xsi, 0);
  LobattoShapes::getShape(n2, dn2, xsi, 1);
  LobattoShapes::getShape(s1, ds1, eta, 0);
  LobattoShapes::getShape(s2, ds2, eta, 1);

  indexType p2 = 1 + order;
  shape(0) = n1;
  shape(1) = n2 * s2;
  shape(p2) = n2 * s1;

  shapeDerivative(0, 0) = dn1;
  shapeDerivative(0, 1) = dn2 * s2;
  shapeDerivative(0, p2) = dn2 * s1;

  shapeDerivative(1, 0) = n1;
  shapeDerivative(1, 1) = n2 * ds2;
  shapeDerivative(1, p2) = n2 * ds1;
  // addition weiterer funktionen
}
void ScaledBoundary2DData::getH1ShapesInternal(
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
      shapeDerivative(0, i + rowcount) = sj * dni;
      shapeDerivative(1, i + rowcount) = ni * dsj;
      LobattoShapes::getShape(sj, dsj, eta, 1);
      shape(i + rowcount + order) = ni * sj;
      shapeDerivative(0, i + rowcount + order) = dni * sj;
      shapeDerivative(1, i + rowcount + order) = ni * dsj;

      // higher funkt
      for (auto j = 0; j < order - 1; ++j) {
        gradeta = j + 2;
        LobattoShapes::getShape(sj, dsj, eta, gradeta);
        counter = 1 + i + j + rowcount;
        shape(counter) = ni * sj;
        shapeDerivative(0, counter) = dni * sj;
        shapeDerivative(1, counter) = ni * dsj;
      }

      rowcount += order;
    }
  } else {
    shape.resize(0);
    shapeDerivative.resize(2, 0);
  }
}

auto ScaledBoundary2DData::getH1Shapes(indexType order,
                                       IntegrationPoint &IntegrationPt)
    -> H1Shapes {

  indexType nRand = order * this->m_edges.size();
  indexType nInnen = nRand * (order - 1);
  indexType numshapes = nRand + nInnen + 1;

  H1Shapes shapes(numshapes, 2);
  shapes.shapes.setZero();
  shapes.shapeDeriv.setZero();

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

  if ((IntegrationPt.sectionNumber + 1) != this->m_edges.size()) {
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
          this->m_edges.size() + IntegrationPt.sectionNumber * (order - 1);
      {
        auto &tempEdge = *m_edges_pointers[IntegrationPt.sectionNumber];

        tempEdge.getH1ShapesInternal(order, tempshape, tempEdgeDerivative, eta);

        for (auto i = 0; i < tempshape.rows(); ++i) {
          shapes.shapes(counter) = tempshape(i) * n2;
          shapes.shapeDeriv(0, counter) = tempshape(i) * dn2;
          shapes.shapeDeriv(1, counter) = tempEdgeDerivative(i) * n2;
          ++counter;
        }
      }
      // Internal
      this->getH1ShapesInternal(order, tempshape, tempshapeDerivative, xsi,
                                eta);
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

    // addition weiterer funktionen
    if (order > 1) {
      Types::VectorX<prec> tempshape;
      Types::VectorX<prec> tempEdgeDerivative;
      Types::Matrix2X<prec> tempshapeDerivative;
      // Edge
      indexType counter =
          this->m_edges.size() + IntegrationPt.sectionNumber * (order - 1);
      {
        auto &tempEdge = *m_edges_pointers[IntegrationPt.sectionNumber];

        tempEdge.getH1ShapesInternal(order, tempshape, tempEdgeDerivative, eta);

        for (auto i = 0; i < tempshape.rows(); ++i) {
          shapes.shapes(counter) = tempshape(i) * n2;
          shapes.shapeDeriv(0, counter) = tempshape(i) * dn2;
          shapes.shapeDeriv(1, counter) = tempEdgeDerivative(i) * n2;
          ++counter;
        }
      }
      // Internal
      this->getH1ShapesInternal(order, tempshape, tempshapeDerivative, xsi,
                                eta);
      counter = nRand + IntegrationPt.sectionNumber * order;
      indexType rowcounter = 0;
      for (auto i = 0; i < tempshape.rows(); ++i) {
        if ((i + 1) % (order + 1) == 0) {
          shapes.shapes(nRand + rowcounter) = tempshape(i);
          shapes.shapeDeriv(0, nRand + rowcounter) = tempshapeDerivative(0, i);
          shapes.shapeDeriv(1, nRand + rowcounter) = tempshapeDerivative(1, i);
          rowcounter += nRand;
          counter += nRand - order;
        } else {
          shapes.shapes(counter) = tempshape(i);
          shapes.shapeDeriv(0, counter) = tempshapeDerivative(0, i);
          shapes.shapeDeriv(1, counter) = tempshapeDerivative(1, i);
          ++counter;
        }
      }
    }
  }

  return shapes;
}

// HDivShapes
void ScaledBoundary2DData::setHDivShapes(indexType meshId, indexType order,
                                         NodeTypes type) {}

void ScaledBoundary2DData::getHDivDofs(std::vector<DegreeOfFreedom *> &Dofs,
                                       indexType meshID, indexType order,
                                       NodeTypes type) {}

void ScaledBoundary2DData::getHDivShapes(indexType order,
                                         Types::Matrix2X<prec> &shape,
                                         Types::VectorX<prec> &dshape, prec xi,
                                         prec eta) {}

// L2 Shapes

void ScaledBoundary2DData::getL2Dofs(std::vector<DegreeOfFreedom *> &Dofs,
                                     indexType meshID, indexType order) {
  auto nodeList = this->getNodeSetNodeListMeshId(meshID);
  nodeList.addDofsToVector(Dofs);
}

void ScaledBoundary2DData::setL2Shapes(indexType meshId, indexType order,
                                       NodeTypes type) {
  if (order == 0) {
    this->setNodeSet(meshId, 1, type);
  }
  // nur Knoten an Kante und im Inneren,Knoten an Vertexe nicht beachtet
  if (order >= 1) {
    // indexType numNodes = 1 + (order - 1) * this->edges.size() + (order - 1) *
    // (this->verts.size() + (order - 1) * this->edges.size());
    indexType numNodes =
        (this->m_verts.size() + (order - 1) * this->m_edges.size()) * order + 1;
    this->setNodeSet(meshId, numNodes, type);
  }
}

auto ScaledBoundary2DData::getL2Shapes(indexType order,
                                       IntegrationPoint &IntegrationPt)
    -> L2Shapes {
  L2Shapes shapes;

  indexType nRand = order * this->m_edges.size();
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

  if ((IntegrationPt.sectionNumber + 1) != this->m_edges.size()) {
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
          this->m_edges.size() + IntegrationPt.sectionNumber * (order - 1);
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
      this->getL2ShapesInternal(order, tempshape, tempshapeDerivative, xsi,
                                eta);
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
          this->m_edges.size() + IntegrationPt.sectionNumber * (order - 1);
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
      this->getL2ShapesInternal(order, tempshape, tempshapeDerivative, xsi,
                                eta);
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

void ScaledBoundary2DData::getL2ShapesInternal(
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

void ScaledBoundary2DData::setAllNodeBoundaryConditionMeshId(indexType meshId,
                                                             indexType dof) {
  GeometryBaseData::setAllNodeBoundaryConditionMeshId(meshId, dof);
  // Fall f�r Druck ausklammern
  if (meshId != 2) {
    for (auto edge : this->m_edges_pointers) {
      edge->setAllNodeBoundaryConditionMeshId(meshId, dof);
    }
  }
}


void ScaledBoundary2DData::flip() {
  throw std::runtime_error(
      "ERROR in ScaledBoundary2D::flip: Method not implemented!");
}

void ScaledBoundary2DData::rotate(indexType n) {
  throw std::runtime_error(
      "ERROR in ScaledBoundary2D::rotate: Method not implemented!");
}

auto ScaledBoundary2DData::computeMeanCoordinate() -> Types::Vector3<prec> {
  throw std::runtime_error("ERROR in ScaledBoundary2D::computeMeanCoordinate: "
                           "Method not implemented!");
  return Types::Vector3<prec>();
}

auto ScaledBoundary2DData::getFaceNormal() -> Types::Vector3<prec> {
  throw std::runtime_error("ERROR in ScaledBoundary2D::getFaceNormal: "
                           "Method not implemented!");
  return Types::Vector3<prec>();
}

void ScaledBoundary2DData::set_geometry_pointers(GeometryData &geoData) {
  m_edges_pointers.clear();
  m_vert_pointers.clear();
  for (auto i : m_verts) {
    m_vert_pointers.push_back(&geoData.getVertexData(i));
  }
  for (auto i : m_edges) {
    m_edges_pointers.push_back(&geoData.getEdgeData(i));
  }
}

auto ScaledBoundary2DData::getVertexNumber(indexType localNumber) -> indexType {
  return m_verts[localNumber];
}

const GeometryTypes ScaledBoundary2DData::type =
    GeometryTypes::ScaledBoundary2D;
} // namespace HierAMuS::Geometry
