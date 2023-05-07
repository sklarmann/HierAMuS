// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <geometry/GeometryData.h>
#include <equations/NodeSet.h>

#include <geometry/Edges.h>
#include <geometry/Faces.h>
#include <geometry/LinearPrism.h>
#include <geometry/Vertex.h>


#include <stdexcept>
#include <types/MatrixTypes.h>

#include <pointercollection/pointercollection.h>
#include "geometry/GeometryData.h"
#include <vector>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <iomanip>

namespace HierAMuS {
namespace Geometry {

LinearPrism::LinearPrism() : Volumes() {}

LinearPrism::~LinearPrism() {}

const GeometryTypes &LinearPrism::getType() { return this->type; }

void LinearPrism::print(PointerCollection &pointers) {
  //
  auto &Logger = pointers.getSPDLogger();

  Logger.debug("Linear Prism id:  {:>12}", this->id);
  Logger.debug("Vertices:  {}", fmt::join(verts, " "));
  Logger.debug("Edges:     {}", fmt::join(edges, " "));
  Logger.debug("Faces:     {}", fmt::join(faces, " "));

  this->printEqInfo(pointers);
}

void LinearPrism::setVerts(GeometryData &geoData, std::vector<indexType> &vertsIn) {
  if (vertsIn.size() == 6) {
    for (auto i = 0; i < 6; ++i) {
      this->verts[i] = vertsIn[i];
    }
  } else {
    throw std::runtime_error("In linear prism geometry element method "
                             "setVerts, wrong amount of vertices passed.");
  }
}

void LinearPrism::getVerts(std::vector<indexType> &vertsOut) {
  vertsOut.assign(std::begin(this->verts), std::end(this->verts));
}

void LinearPrism::getVerts(PointerCollection &pointers,
                           std::vector<Base *> &vertsOut) {
  vertsOut.clear();
  for (auto i = 0; i < 6; ++i) {
    vertsOut.push_back(&pointers.getGeometryData()->getVertex(this->verts[i]));
  }
}

void LinearPrism::setEdges(const std::vector<indexType> &edgesIn) {
  if (edgesIn.size() == 9) {
    for (auto i = 0; i < 9; ++i) {
      this->edges[i] = edgesIn[i];
    }
  } else {
    throw std::runtime_error("In linear prism geometry element method "
                             "setEdges, wrong amount of edges passed.");
  }
}

void LinearPrism::getEdges(std::vector<indexType> &edgesOut) {
  edgesOut.assign(std::begin(this->edges), std::end(this->edges));
}

void LinearPrism::getEdges(PointerCollection &pointers,
                           std::vector<Base *> &edgesOut) {
  edgesOut.clear();
  for (auto i = 0; i < 9; ++i) {
    edgesOut.push_back(&pointers.getGeometryData()->getEdge(this->edges[i]));
  }
}

void LinearPrism::setFaces(std::vector<indexType> &facesIn) {
  if (facesIn.size() == 5) {
    for (auto i = 0; i < 5; ++i) {
      this->faces[i] = facesIn[i];
    }
  } else {
    throw std::runtime_error("In linear prism geometry element method "
                             "setFaces, wrong amount of faces passed.");
  }
}

void LinearPrism::getFaces(std::vector<indexType> &facesOut) {
  facesOut.assign(std::begin(this->faces), std::end(this->faces));
}

void LinearPrism::getJacobian(PointerCollection &pointers,
                              Types::Matrix33<prec> &jacobi, prec xsi,
                              prec eta, prec zeta) {

  Types::Vector3<prec> coors;
  jacobi.setZero();
  Types::Matrix3X<prec> dshape;
  Types::VectorX<prec> shape;
  shape.resize(6);
  dshape.resize(3, 6);
  this->getH1Shapes(pointers, 1, shape, dshape, xsi, eta, zeta);
  for (auto i = 0; i < 6; i++) {
    auto &tempVert = pointers.getGeometryData()->getVertex(this->verts[i]);
    coors = tempVert.getCoordinates();
    for (auto j = 0; j < 3; j++) {
      for (auto k = 0; k < 3; k++) {
        jacobi(k, j) += dshape(j, i) * coors[k];
      }
    }
  }
}

void LinearPrism::setH1Shapes(PointerCollection &pointers,
                              indexType meshId, indexType order,
                              NodeTypes type) {
  Faces *tempFace;
  for (auto i = 0; i < 5; i++) {
    tempFace = pointers.getGeometryData()->getFace(this->faces[i]);
    tempFace->setH1Shapes(pointers, meshId, order, type);
  }
  if (order > 1) {
    this->setH1ShapesInternal(pointers, meshId, order, type);
  }
}

void LinearPrism::setH1ShapesInternal(PointerCollection &pointers,
                                      indexType meshId,
                                      indexType order,
                                      NodeTypes type) {
  if (order > 1) {
    indexType num = order - 1;
    num *= num * num;
    this->setNodeSet(pointers, meshId, num, type);
  }
}

void LinearPrism::getH1Dofs(PointerCollection &pointers,
                            std::vector<DegreeOfFreedom *> &Dofs,
                            indexType meshID, indexType order) {

  NodeSet *tempSet;
  std::vector<DegreeOfFreedom *> tdofs;
  for (auto i = 0; i < 6; i++) {
    auto &tempGeo = pointers.getGeometryData()->getVertex(this->verts[i]);
    tempSet = tempGeo.getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }

  if (order > 1) {
    for (auto i = 0; i < 9; i++) {
      auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[i]);
      tempEdge.getH1DofsInternal(pointers, Dofs, meshID, order);
    }
    Faces *tempFace;
    for (auto i = 0; i < 5; i++) {
      tempFace = pointers.getGeometryData()->getFace(this->faces[i]);
      tempFace->getH1DofsInternal(pointers, Dofs, meshID, order);
    }
    this->getH1DofsInternal(pointers, Dofs, meshID, order);
  }
}

void LinearPrism::getH1DofsInternal(PointerCollection &pointers,
                                    std::vector<DegreeOfFreedom *> &Dofs,
                                    indexType meshID,
                                    indexType order) {
  if (order > 1) {
    NodeSet *tempSet;
    std::vector<DegreeOfFreedom *> tdofs;
    tempSet = this->getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
}

void LinearPrism::getH1Shapes(PointerCollection &pointers,
                              indexType order,
                              Types::VectorX<prec> &shape,
                              Types::Matrix3X<prec> &shapeDerivative,
                              prec xsi, prec eta,
                              prec zeta) {

  throw std::runtime_error(
      "Error in linear prism, method getH1Shapes not implemented!");

  indexType numshapes = (order + 1);
  numshapes *= numshapes * numshapes;

  shape.resize(numshapes);
  shapeDerivative.resize(3, numshapes);


  auto &lshape = pointers.getLobattoShapes();

  auto [s1, ds1] = lshape.getShape(xsi, 0);
  auto [s2, ds2] = lshape.getShape(xsi, 1);

  auto [s3, ds3] = lshape.getShape(eta, 0);
  auto [s4, ds4] = lshape.getShape(eta, 1);

  auto [s5, ds5] = lshape.getShape(zeta, 0);
  auto [s6, ds6] = lshape.getShape(zeta, 1);



  shape(0) = s1 * s3 * s5;
  shapeDerivative(0, 0) = ds1 * s3 * s5;
  shapeDerivative(1, 0) = s1 * ds3 * s5;
  shapeDerivative(2, 0) = s1 * s3 * ds5;

  shape(1) = s2 * s3 * s5;
  shapeDerivative(0, 1) = ds2 * s3 * s5;
  shapeDerivative(1, 1) = s2 * ds3 * s5;
  shapeDerivative(2, 1) = s2 * s3 * ds5;

  shape(2) = s2 * s4 * s5;
  shapeDerivative(0, 2) = ds2 * s4 * s5;
  shapeDerivative(1, 2) = s2 * ds4 * s5;
  shapeDerivative(2, 2) = s2 * s4 * ds5;

  shape(3) = s1 * s4 * s5;
  shapeDerivative(0, 3) = ds1 * s4 * s5;
  shapeDerivative(1, 3) = s1 * ds4 * s5;
  shapeDerivative(2, 3) = s1 * s4 * ds5;

  shape(4) = s1 * s3 * s6;
  shapeDerivative(0, 4) = ds1 * s3 * s6;
  shapeDerivative(1, 4) = s1 * ds3 * s6;
  shapeDerivative(2, 4) = s1 * s3 * ds6;

  shape(5) = s2 * s3 * s6;
  shapeDerivative(0, 5) = ds2 * s3 * s6;
  shapeDerivative(1, 5) = s2 * ds3 * s6;
  shapeDerivative(2, 5) = s2 * s3 * ds6;

  shape(6) = s2 * s4 * s6;
  shapeDerivative(0, 6) = ds2 * s4 * s6;
  shapeDerivative(1, 6) = s2 * ds4 * s6;
  shapeDerivative(2, 6) = s2 * s4 * ds6;

  shape(7) = s1 * s4 * s6;
  shapeDerivative(0, 7) = ds1 * s4 * s6;
  shapeDerivative(1, 7) = s1 * ds4 * s6;
  shapeDerivative(2, 7) = s1 * s4 * ds6;

  if (order > 1) {

    Types::VectorX<prec> tempshape;
    Types::Matrix3X<prec> tempdshape;
    this->getH1ShapesInternal(pointers, order, tempshape, tempdshape, xsi, eta,
                              zeta);
  }
}

void LinearPrism::getH1ShapesInternal(PointerCollection &pointers,
                                      indexType order,
                                      Types::VectorX<prec> &shape,
                                      Types::Matrix3X<prec> &shapeDerivative,
                                      prec xsi, prec eta,
                                      prec zeta) {
  if (order > 1) {
    throw std::runtime_error(
        "Higher order shapes not implemented for linear Prism element!");
  }
}

void LinearPrism::getFaces(PointerCollection &pointers,
                           std::vector<Base *> &facesOut) {
  facesOut.clear();
  for (auto i = 0; i < 5; ++i) {
    facesOut.push_back(pointers.getGeometryData()->getFace(this->faces[i]));
  }
}

const GeometryTypes LinearPrism::type = GeometryTypes::LinearPrism;
} // namespace Geometry
} // namespace HierAMuS
