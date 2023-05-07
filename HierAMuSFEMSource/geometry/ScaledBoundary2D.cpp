// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "MatrixTypes.h"
#include "geometry/GeometryTypes.h"

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <geometry/ScaledBoundary2D.h>
#include <pointercollection/pointercollection.h>
#include <sstream>

#include <equations/NodeSet.h>
#include <geometry/Edges.h>
#include <geometry/GeometryData.h>
#include <geometry/Vertex.h>
#include "geometry/GeometryData.h"

#include <vtkCellType.h>

#include "HelperFunctions.h"

namespace HierAMuS::Geometry {

ScaledBoundary2D::ScaledBoundary2D() = default;

ScaledBoundary2D::~ScaledBoundary2D() = default;

auto ScaledBoundary2D::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::ScaledBoundary2D::type;
}

void ScaledBoundary2D::setVerts(GeometryData &geoData, std::vector<indexType> &vertsIn) {
  indexType numVerts = vertsIn.size();
  if (numVerts >= 3) {
    this->verts.resize(numVerts);
    for (auto i = 0; i < vertsIn.size(); ++i) {
      this->verts[i] = vertsIn[i];
    }
  } else {
    std::stringstream temp;
    temp
        << "ScaledBoundary2D geometry element needs atleast 3 vertices. Cannot "
        << "add the given amount of " << vertsIn.size() << " vertices!";
    std::runtime_error(temp.str());
  }
}

void ScaledBoundary2D::setEdges(const std::vector<indexType> &edgesIn) {
  indexType numEdges = edgesIn.size();
  if (edgesIn.size() >= 3) {
    this->edges.resize(numEdges);
    for (auto i = 0; i < edgesIn.size(); ++i) {
      this->edges[i] = edgesIn[i];
    }
  } else {
    std::stringstream temp;
    temp
        << "ScaledBoundary2D geometry element needs atleast 3 vertices. Cannot "
        << "add the given amount of " << edgesIn.size() << " edges!";
    std::runtime_error(temp.str());
  }
}

void ScaledBoundary2D::getVerts(std::vector<indexType> &vertsOut) {
  vertsOut.resize(this->edges.size());
  for (auto i = 0; i < this->edges.size(); ++i) {
    vertsOut[i] = this->verts[i];
  }
}

void ScaledBoundary2D::getVerts(PointerCollection &pointers,
                                std::vector<Base *> &vertsOut) {
  vertsOut.clear();
  for (auto i = 0; i < this->edges.size(); i++) {
    vertsOut.push_back(&pointers.getGeometryData()->getVertex(this->verts[i]));
  }
}

void ScaledBoundary2D::getEdges(std::vector<indexType> &edgesOut) {
  edgesOut.resize(this->edges.size());
  for (auto i = 0; i < this->edges.size(); ++i) {
    edgesOut[i] = this->edges[i];
  }
}
void ScaledBoundary2D::setScalingCenter(prec x, prec y) {
  this->scalingcoor.resize(2);
  this->scalingcoor[0] = x;
  this->scalingcoor[1] = y;
}

void ScaledBoundary2D::computeScalingCenter(PointerCollection &pointers) {
  // this->scalingCenter =
  // pointers.getGeometryData()->requestNewGeometryObject(GeometryTypes::Vertex);
  // auto &vertex = pointers.getGeometryData()->getVertex(this->scalingCenter);
  // vertex.setCoordinates(1,2,0);

  // this->scalingcoor.resize(2);
  // this->scalingcoor[0] = 10;
  // this->scalingcoor[1] = 10;

  // ScalingCenter im Kreuz f�r Viereck

  auto tempcoor0 =
      pointers.getGeometryData()->getVertex(this->verts[0]).getCoordinates();
  auto tempcoor1 =
      pointers.getGeometryData()->getVertex(this->verts[1]).getCoordinates();
  auto tempcoor2 =
      pointers.getGeometryData()->getVertex(this->verts[2]).getCoordinates();
  auto tempcoor3 =
      pointers.getGeometryData()->getVertex(this->verts[3]).getCoordinates();

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

void ScaledBoundary2D::getEdges(PointerCollection &pointers,
                                std::vector<Base *> &edgesOut) {
  edgesOut.clear();
  for (auto edge : this->edges) {
    edgesOut.push_back(&pointers.getGeometryData()->getEdge(edge));
  }
}

auto ScaledBoundary2D::getEdge(indexType local_number) -> indexType {
  return edges[local_number];
}

auto ScaledBoundary2D::hasVertices(indexType v1, indexType v2, indexType v3)
    -> bool {
  if (contains(this->verts, v1)) {
    if (contains(this->verts, v2)) {
      if (contains(this->verts, v3)) {
        return true;
      }
    }
  }
  return false;
}

inline void ScaledBoundary2D::print(PointerCollection &pointers) {
  auto &Logger = pointers.getSPDLogger();
  Logger.debug("Scaled Boundary 2D Face id:   {:>10}",this->id);
  Logger.debug("Vertices:  {}",fmt::join(this->verts," "));
  Logger.debug("Edges:     {}",fmt::join(this->edges," "));
  this->printEqInfo(pointers);
}

auto ScaledBoundary2D::getCoordinates(PointerCollection &pointers, prec xi,
                                      prec eta) -> Types::Vector3<prec> {
  return {};
}

auto ScaledBoundary2D::getIntegrationPoints(PointerCollection &pointers, indexType elementId)
-> IntegrationPoints {
  IntegrationPoints temp = HierAMuS::PointerCollection::getIntegrationPoints(elementId);
  temp.setType(IntegrationType::Scaled2D);
  temp.setNumberOfSections(this->edges.size());
  return temp;
}

auto ScaledBoundary2D::getJacobian(PointerCollection &pointers,
                                   IntegrationPoint &IntegrationPt)
    -> Types::MatrixXX<prec> {

  Types::MatrixXX<prec> jacobi;
  jacobi.resize(2, 2);
  jacobi.setZero();

  auto shapes = this->getH1Shapes(pointers, 1, IntegrationPt);
  indexType currsection = IntegrationPt.sectionNumber;
  int numshapes = shapes.shapes.rows();
  jacobi(0, 0) += shapes.shapeDeriv(0, numshapes - 1) * this->scalingcoor[0];
  jacobi(0, 1) += shapes.shapeDeriv(1, numshapes - 1) * this->scalingcoor[0];
  jacobi(1, 0) += shapes.shapeDeriv(0, numshapes - 1) * this->scalingcoor[1];
  jacobi(1, 1) += shapes.shapeDeriv(1, numshapes - 1) * this->scalingcoor[1];

  if ((currsection + 1) == edges.size()) {
    // int a = currsection + 1;
    auto coord = pointers.getGeometryData()
                     ->getVertex(this->verts[currsection])
                     .getCoordinates();
    jacobi(0, 0) += shapes.shapeDeriv(0, currsection) * coord(0);
    jacobi(0, 1) += shapes.shapeDeriv(1, currsection) * coord(0);
    jacobi(1, 0) += shapes.shapeDeriv(0, currsection) * coord(1);
    jacobi(1, 1) += shapes.shapeDeriv(1, currsection) * coord(1);

    coord =
        pointers.getGeometryData()->getVertex(this->verts[0]).getCoordinates();
    jacobi(0, 0) += shapes.shapeDeriv(0, 0) * coord(0);
    jacobi(0, 1) += shapes.shapeDeriv(1, 0) * coord(0);
    jacobi(1, 0) += shapes.shapeDeriv(0, 0) * coord(1);
    jacobi(1, 1) += shapes.shapeDeriv(1, 0) * coord(1);
  } else {
    for (auto i = currsection; i < currsection + 2; ++i) {

      auto coord = pointers.getGeometryData()
                       ->getVertex(this->verts[i])
                       .getCoordinates();

      jacobi(0, 0) += shapes.shapeDeriv(0, i) * coord(0);
      jacobi(0, 1) += shapes.shapeDeriv(1, i) * coord(0);
      jacobi(1, 0) += shapes.shapeDeriv(0, i) * coord(1);
      jacobi(1, 1) += shapes.shapeDeriv(1, i) * coord(1);
    }
  }
  return jacobi;
}

// H1Shapes

void ScaledBoundary2D::setH1Shapes(PointerCollection &pointers,
                                   indexType meshId, indexType order,
                                   NodeTypes type) {
  for (auto edge : this->edges) {
    auto &edgeTemp = pointers.getGeometryData()->getEdge(edge);
    edgeTemp.setH1Shapes(pointers, meshId, order, type);
  }
  this->setH1ShapesInternal(pointers, meshId, order, type);
}

void ScaledBoundary2D::setH1ShapesInternal(PointerCollection &pointers,
                                           indexType meshId, indexType order,
                                           NodeTypes type) {
  // ScalingCenter
  if (order == 1) {
    this->setNodeSet(pointers, meshId, 1, type);
  }

  if (order > 1) {
    indexType numNodes = 1 + (order - 1) * (this->verts.size() +
                                            (order - 1) * this->edges.size());

    this->setNodeSet(pointers, meshId, numNodes, type);
  }
}

void ScaledBoundary2D::getH1Dofs(PointerCollection &pointers,
                                 std::vector<DegreeOfFreedom *> &Dofs,
                                 indexType meshID, indexType order) {
  std::vector<DegreeOfFreedom *> tdofs;
  NodeSet *tempSet;
  for (auto i = 0; i < this->edges.size(); ++i) {
    auto &tempVert = pointers.getGeometryData()->getVertex(this->verts[i]);
    tempSet = tempVert.getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
  this->getH1DofsInternal(pointers, Dofs, meshID, order);
}

void ScaledBoundary2D::getH1DofsInternal(PointerCollection &pointers,
                                         std::vector<DegreeOfFreedom *> &Dofs,
                                         indexType meshID, indexType order) {
  // ScalingCenter
  if (order == 1) {
    std::vector<DegreeOfFreedom *> tdofs;
    NodeSet *tempSet;
    tempSet = this->getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
  if (order > 1) {
    std::vector<DegreeOfFreedom *> tdofs;
    NodeSet *tempSet;
    for (auto edge : this->edges) {
      auto &tempEdge = pointers.getGeometryData()->getEdge(edge);
      tempEdge.getH1DofsInternal(pointers, Dofs, meshID, order);
    }
    // innere Nodes
    tempSet = this->getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
}

void ScaledBoundary2D::getH1Shapes(PointerCollection &pointers, indexType order,
                                   Types::VectorX<prec> &shape,
                                   Types::Matrix2X<prec> &shapeDerivative,
                                   prec xsi, prec eta) {
  // fehlerbehaftet

  indexType nRand = order * this->edges.size();
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
void ScaledBoundary2D::getH1ShapesInternal(
    PointerCollection &pointers, indexType order, Types::VectorX<prec> &shape,
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

auto ScaledBoundary2D::getH1Shapes(PointerCollection &pointers, indexType order,
                                   IntegrationPoint &IntegrationPt)
    -> H1Shapes {
  H1Shapes shapes;

  indexType nRand = order * this->edges.size();
  indexType nInnen = nRand * (order - 1);
  indexType numshapes = nRand + nInnen + 1;

  shapes.shapes.resize(numshapes);
  shapes.shapes.setZero();
  shapes.shapeDeriv.resize(2, numshapes);
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
        auto &tempEdge = pointers.getGeometryData()->getEdge(
            this->edges[IntegrationPt.sectionNumber]);

        tempEdge.getH1ShapesInternal(pointers, order, tempshape,
                                     tempEdgeDerivative, eta);

        for (auto i = 0; i < tempshape.rows(); ++i) {
          shapes.shapes(counter) = tempshape(i) * n2;
          shapes.shapeDeriv(0, counter) = tempshape(i) * dn2;
          shapes.shapeDeriv(1, counter) = tempEdgeDerivative(i) * n2;
          ++counter;
        }
      }
      // Internal
      this->getH1ShapesInternal(pointers, order, tempshape, tempshapeDerivative,
                                xsi, eta);
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
          this->edges.size() + IntegrationPt.sectionNumber * (order - 1);
      {
        auto &tempEdge = pointers.getGeometryData()->getEdge(
            this->edges[IntegrationPt.sectionNumber]);

        tempEdge.getH1ShapesInternal(pointers, order, tempshape,
                                     tempEdgeDerivative, eta);

        for (auto i = 0; i < tempshape.rows(); ++i) {
          shapes.shapes(counter) = tempshape(i) * n2;
          shapes.shapeDeriv(0, counter) = tempshape(i) * dn2;
          shapes.shapeDeriv(1, counter) = tempEdgeDerivative(i) * n2;
          ++counter;
        }
      }
      // Internal
      this->getH1ShapesInternal(pointers, order, tempshape, tempshapeDerivative,
                                xsi, eta);
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
void ScaledBoundary2D::setHDivShapes(PointerCollection &pointers,
                                     indexType meshId, indexType order,
                                     NodeTypes type) {}

void ScaledBoundary2D::getHDivDofs(PointerCollection &pointers,
                                   std::vector<DegreeOfFreedom *> &Dofs,
                                   indexType meshID, indexType order,
                                   NodeTypes type) {}

void ScaledBoundary2D::getHDivShapes(PointerCollection &pointers,
                                     indexType order,
                                     Types::Matrix2X<prec> &shape,
                                     Types::VectorX<prec> &dshape,
                                     prec xi, prec eta) {}

// L2 Shapes

void ScaledBoundary2D::getL2Dofs(PointerCollection &pointers,
                                 std::vector<DegreeOfFreedom *> &Dofs,
                                 indexType meshID, indexType order) {
  std::vector<DegreeOfFreedom *> tdofs;
  NodeSet *tempSet;
  tempSet = this->getSetMeshId(pointers, meshID);
  tdofs = tempSet->getDegreesOfFreedom(pointers);
  Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
}

void ScaledBoundary2D::setL2Shapes(PointerCollection &pointers,
                                   indexType meshId, indexType order,
                                   NodeTypes type) {
  if (order == 0) {
    this->setNodeSet(pointers, meshId, 1, type);
  }
  // nur Knoten an Kante und im Inneren,Knoten an Vertexe nicht beachtet
  if (order >= 1) {
    // indexType numNodes = 1 + (order - 1) * this->edges.size() + (order - 1) *
    // (this->verts.size() + (order - 1) * this->edges.size());
    indexType numNodes =
        (this->verts.size() + (order - 1) * this->edges.size()) * order + 1;
    this->setNodeSet(pointers, meshId, numNodes, type);
  }
}

auto ScaledBoundary2D::getL2Shapes(PointerCollection &pointers, indexType order,
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
      this->getL2ShapesInternal(pointers, order, tempshape, tempshapeDerivative,
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
      this->getL2ShapesInternal(pointers, order, tempshape, tempshapeDerivative,
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

void ScaledBoundary2D::getL2ShapesInternal(
    PointerCollection &pointers, indexType order, Types::VectorX<prec> &shape,
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

void ScaledBoundary2D::setAllNodeBoundaryConditionMeshId(
    PointerCollection &pointers, indexType meshId, indexType dof) {
  Base::setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);
  // Fall f�r Druck ausklammern
  if (meshId != 2) {
    for (auto edge : this->edges) {
      auto &tempEdge = pointers.getGeometryData()->getEdge(edge);
      tempEdge.setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);
    }
  }
}

void ScaledBoundary2D::geometryToParaview(PointerCollection &pointers,
                                          vtkPlotInterface &paraviewAdapter,
                                          indexType mainMesh,
                                          indexType subMesh) {
  indexType numPoints = this->edges.size();
  std::vector<indexType> points(numPoints);
  points.clear();
  for (auto i : this->verts) {
    auto &vert = pointers.getGeometryData()->getVertex(i);
    vert.geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
    points.push_back(i);
  }

  paraviewAdapter.addCell(mainMesh, subMesh, this->id, 1, points, numPoints,
                          VTK_POLYGON);
}

void ScaledBoundary2D::computeWeightsParaview(PointerCollection &pointers,
                                              vtkPlotInterface &paraviewAdapter,
                                              indexType mainMesh,
                                              indexType subMesh) {

  auto GP = this->getIntegrationPoints(pointers,-1);
  GP.setOrder(2);

  for (auto i : GP) {
    auto jaco = this->getJacobian(pointers, i);
    auto shapes = this->getH1Shapes(pointers, 1, i);
    prec dA = jaco.determinant() * i.weight;

    for (auto i = 0; i < 4; ++i) {
      auto &vert = pointers.getGeometryData()->getVertex(this->verts[i]);
      std::vector<prec> val;
      val.push_back(shapes.shapes(i) * dA);
      paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, val, vert.getId(),
                                           1, paraviewNames::weightName());
    }
  }
}

void ScaledBoundary2D::H1SolutionToParaview(PointerCollection &pointers,
                                            vtkPlotInterface &paraviewAdapter,
                                            indexType mainMesh,
                                            indexType subMesh, indexType meshId,
                                            indexType order,
                                            std::string &name) {
  std::vector<DegreeOfFreedom *> Dofs;
  this->getH1Dofs(pointers, Dofs, meshId, order);
  auto solution = pointers.getSolutionState()->getSolution(Dofs);
  for (auto i = 0; i < this->edges.size(); ++i) {
    auto &V = pointers.getGeometryData()->getVertex(this->verts[i]);
    std::vector<prec> sol(3);
    for (auto j = 0; j < 3; ++j) {
      sol[j] = solution(3 * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, V.getId(), sol, 3, name);
  }
}
void ScaledBoundary2D::H1DataToParaview(PointerCollection &pointers,
                                        vtkPlotInterface &paraviewAdapter,
                                        indexType mainMesh, indexType subMesh,
                                        Types::VectorX<prec> &Data,
                                        indexType numberComponents,
                                        indexType order, std::string &name) {
  for (auto i = 0; i < this->edges.size(); ++i) {
    auto &V = pointers.getGeometryData()->getVertex(this->verts[i]);
    std::vector<prec> sol(numberComponents);
    for (auto j = 0; j < numberComponents; ++j) {
      sol[j] = Data(numberComponents * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, V.getId(), sol,
                                 numberComponents, name);
  }
}

void ScaledBoundary2D::projectDataToParaviewVertices(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {

  auto shapes = this->getH1Shapes(pointers, 1, IntegrationPt);
  std::vector<prec> vals(numberComponents);

  auto jaco = this->getJacobian(pointers, IntegrationPt);
  auto dA = jaco.determinant() * IntegrationPt.weight;
  for (auto i = 0; i < this->edges.size(); ++i) {
    auto &V = pointers.getGeometryData()->getVertex(this->verts[i]);
    for (auto j = 0; j < numberComponents; ++j) {
      vals[j] = data(j) * shapes.shapes(i) * dA;
    }
    paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, vals, V.getId(),
                                         numberComponents, name);
  }
}

void ScaledBoundary2D::flip() {
  throw std::runtime_error(
      "ERROR in ScaledBoundary2D::flip: Method not implemented!");
}

void ScaledBoundary2D::rotate(indexType n) {
  throw std::runtime_error(
      "ERROR in ScaledBoundary2D::rotate: Method not implemented!");
}

auto ScaledBoundary2D::computeMeanCoordinate(PointerCollection &pointers)
    -> Types::Vector3<prec> {
  throw std::runtime_error(
      "ERROR in ScaledBoundary2D::computeMeanCoordinate: Method not implemented!");
  return Types::Vector3<prec>();
}

const GeometryTypes ScaledBoundary2D::type = GeometryTypes::ScaledBoundary2D;
} // namespace HierAMuS::Geometry
