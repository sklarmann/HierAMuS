// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "geometry/Base.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

#include <equations/NodeSet.h>
#include <geometry/GeometryData.h>

#include <geometry/Edges.h>
#include <geometry/Faces.h>
#include <geometry/LinearBrick.h>
#include <geometry/Vertex.h>

#include <types/MatrixTypes.h>

#include "geometry/GeometryData.h"
#include <pointercollection/pointercollection.h>
#include <vector>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <iomanip>

#include "plot/vtkplotClass.h"

#include <vtkCellType.h>
namespace HierAMuS::Geometry {

LinearBrick::LinearBrick()
    : m_edges({-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}),
      m_faces({-1, -1, -1, -1, -1, -1}) {}

LinearBrick::~LinearBrick() = default;

auto LinearBrick::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::LinearBrick::type;
}

void LinearBrick::print(PointerCollection &pointers) {
  //
  auto &Logger = pointers.getSPDLogger();

  Logger.debug("Linear Brick id: {:8d}", this->id);
  Logger.debug("Vertices: {}", fmt::join(m_verts, " "));
  Logger.debug("Edges:    {}", fmt::join(m_edges, " "));
  Logger.debug("Faces:    {}", fmt::join(m_faces, " "));

  this->printEqInfo(pointers);
}

void LinearBrick::setVerts(GeometryData &geoData,
                           std::vector<indexType> &vertsIn) {
  if (vertsIn.size() == 8) {
    for (auto i = 0; i < 8; ++i) {
      this->m_verts[i] = vertsIn[i];
    }
  } else {
    // TODO throw exception
  }
}

void LinearBrick::setAllNodeBoundaryConditionMeshId(PointerCollection &pointers,
                                                    indexType meshId,
                                                    indexType dof) {
  Base::setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);

  for (auto i : m_faces) {
    auto tempFace = pointers.getGeometryData()->getFace(i);
    tempFace->setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);
  }
}

void LinearBrick::getVerts(std::vector<indexType> &vertsOut) {
  vertsOut.clear();
  vertsOut.insert(vertsOut.end(), m_verts.begin(), m_verts.end());
  // vertsOut.assign(this->m_verts, this->m_verts + 8);
}

void LinearBrick::getVerts(PointerCollection &pointers,
                           std::vector<Base *> &vertsOut) {
  vertsOut.clear();
  for (auto i = 0; i < 8; ++i) {
    vertsOut.push_back(
        &pointers.getGeometryData()->getVertex(this->m_verts[i]));
  }
}

auto LinearBrick::getCoordinates(PointerCollection &pointers,
                                 IntegrationPoint &IntPoint)
    -> Types::Vector3<prec> {

  auto shp = this->getH1Shapes(pointers, 1, IntPoint);

  Types::Vector3<prec> coor;
  coor.setZero();

  for (indexType i = 0; i < 8; ++i) {
    auto &Vert = pointers.getGeometryData()->getVertex(this->m_verts[i]);
    coor += shp.shapes(i) * Vert.getCoordinates(pointers, IntPoint);
  }

  return coor;
}

void LinearBrick::setEdges(const std::vector<indexType> &edgesIn) {
  if (edgesIn.size() == 12) {
    for (auto i = 0; i < 12; ++i) {
      this->m_edges[i] = edgesIn[i];
    }
  } else {
    // TODO throw exception
  }
}

void LinearBrick::getEdges(std::vector<indexType> &edgesOut) {
  edgesOut.clear();
  edgesOut.insert(edgesOut.end(), m_edges.begin(), m_edges.end());
  // edgesOut.assign(this->m_edges, this->m_edges + 12);
}

void LinearBrick::getEdges(PointerCollection &pointers,
                           std::vector<Base *> &edgesOut) {
  edgesOut.clear();
  for (auto i = 0; i < 12; ++i) {
    edgesOut.push_back(&pointers.getGeometryData()->getEdge(this->m_edges[i]));
  }
}

void LinearBrick::setFaces(std::vector<indexType> &facesIn) {
  if (facesIn.size() == 6) {
    for (auto i = 0; i < 6; ++i) {
      this->m_faces[i] = facesIn[i];
    }
  } else {
    // TODO throw exception
  }
}

void LinearBrick::getFaces(std::vector<indexType> &facesOut) {
  facesOut.clear();
  facesOut.insert(facesOut.end(), m_faces.begin(), m_faces.end());
  // facesOut.assign(this->m_faces, this->m_faces + 6);
}

auto LinearBrick::getIntegrationPoints(PointerCollection &pointers,
                                       indexType elementId)
    -> IntegrationPoints {
  auto points = HierAMuS::PointerCollection::getIntegrationPoints(elementId);
  points.setType(IntegrationType::Gauss3D);
  return points;
}

void LinearBrick::getJacobian(PointerCollection &pointers,
                              Types::Matrix33<prec> &jacobi, prec xsi, prec eta,
                              prec zeta) {

  Types::Vector3<prec> coors;
  jacobi.setZero();
  Types::Matrix3X<prec> dshape;
  Types::VectorX<prec> shape;
  shape.resize(8);
  dshape.resize(3, 8);
  this->getH1Shapes(pointers, 1, shape, dshape, xsi, eta, zeta);
  auto geoData = pointers.getGeometryData();
  for (auto i = 0; i < 8; i++) {
    auto coors = geoData->getVertex(i).getCoordinates();
    for (auto j = 0; j < 3; j++) {
      for (auto k = 0; k < 3; k++) {
        jacobi(k, j) += dshape(j, i) * coors[k];
      }
    }
  }
}

auto LinearBrick::getJacobian(PointerCollection &pointers,
                              IntegrationPoint &IntegrationPt)
    -> Types::MatrixXX<prec> {
  Types::MatrixXX<prec> jacobi(3, 3);
  jacobi.setZero();
  auto shapes = this->getH1Shapes(pointers, 1, IntegrationPt);
  indexType cc = 0;
  auto geoData = pointers.getGeometryData();
  for (auto i : this->m_verts) {
    auto coors = geoData->getVertex(i).getCoordinates();
    for (auto j = 0; j < 3; j++) {
      for (auto k = 0; k < 3; k++) {
        jacobi(k, j) += shapes.shapeDeriv(j, cc) * coors(k);
      }
    }
    ++cc;
  }
  return jacobi;
}

void LinearBrick::setH1Shapes(PointerCollection &pointers, indexType meshId,
                              indexType order, NodeTypes type) {
  Faces *tempFace;
  for (auto i = 0; i < 6; i++) {
    tempFace = pointers.getGeometryData()->getFace(this->m_faces[i]);
    tempFace->setH1Shapes(pointers, meshId, order, type);
  }
  if (order > 1) {
    this->setH1ShapesInternal(pointers, meshId, order, type);
  }
}

void LinearBrick::setH1ShapesInternal(PointerCollection &pointers,
                                      indexType meshId, indexType order,
                                      NodeTypes type) {
  if (order > 1) {
    indexType num = order - 1;
    num *= num * num;
    this->setNodeSet(pointers, meshId, num, type);
  }
}

void LinearBrick::getH1Dofs(PointerCollection &pointers,
                            std::vector<DegreeOfFreedom *> &Dofs,
                            indexType meshID, indexType order) {

  NodeSet *tempSet;
  std::vector<DegreeOfFreedom *> tdofs;
  auto geoData = pointers.getGeometryData();
  for (auto i = 0; i < 8; i++) {
    auto &tempGeo = geoData->getVertex(this->m_verts[i]);
    tempSet = tempGeo.getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }

  if (order > 1) {
    for (auto i = 0; i < 12; i++) {
      auto &tempEdge = geoData->getEdge(this->m_edges[i]);
      tempEdge.getH1DofsInternal(pointers, Dofs, meshID, order);
    }

    for (auto i = 0; i < 6; i++) {
      auto tempFace = geoData->getFace(this->m_faces[i]);
      tempFace->getH1DofsInternal(pointers, Dofs, meshID, order);
    }
    this->getH1DofsInternal(pointers, Dofs, meshID, order);
  }
}

void LinearBrick::getH1DofsInternal(PointerCollection &pointers,
                                    std::vector<DegreeOfFreedom *> &Dofs,
                                    indexType meshID, indexType order) {
  if (order > 1) {
    NodeSet *tempSet;
    std::vector<DegreeOfFreedom *> tdofs;
    tempSet = this->getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
}

void LinearBrick::getH1Shapes(PointerCollection &pointers, indexType order,
                              Types::VectorX<prec> &shape,
                              Types::Matrix3X<prec> &shapeDerivative, prec xsi,
                              prec eta, prec zeta) {
  indexType numshapes = (order + 1);
  numshapes *= numshapes * numshapes;

  shape.resize(numshapes);
  shapeDerivative.resize(3, numshapes);

  auto [s1, ds1] = HierAMuS::LobattoShapes::getShape(xsi, 0);
  auto [s2, ds2] = HierAMuS::LobattoShapes::getShape(xsi, 1);

  auto [s3, ds3] = HierAMuS::LobattoShapes::getShape(eta, 0);
  auto [s4, ds4] = HierAMuS::LobattoShapes::getShape(eta, 1);

  auto [s5, ds5] = HierAMuS::LobattoShapes::getShape(zeta, 0);
  auto [s6, ds6] = HierAMuS::LobattoShapes::getShape(zeta, 1);

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

void LinearBrick::getH1ShapesInternal(PointerCollection &pointers,
                                      indexType order,
                                      Types::VectorX<prec> &shape,
                                      Types::Matrix3X<prec> &shapeDerivative,
                                      prec xsi, prec eta, prec zeta) {
  if (order > 1) {
    throw std::runtime_error(
        "Higher order shapes not implemented for Brick element!");
  }
}

auto LinearBrick::getH1Shapes(PointerCollection &pointers, indexType order,
                              IntegrationPoint &IntegrationPt) -> H1Shapes {
  H1Shapes shapes;

  indexType numshapes = (order + 1);
  numshapes *= numshapes * numshapes;

  shapes.shapes.resize(numshapes);
  shapes.shapeDeriv.resize(3, numshapes);

  auto [s1, ds1] = HierAMuS::LobattoShapes::getShape(IntegrationPt.xi, 0);
  auto [s2, ds2] = HierAMuS::LobattoShapes::getShape(IntegrationPt.xi, 1);

  auto [s3, ds3] = HierAMuS::LobattoShapes::getShape(IntegrationPt.eta, 0);
  auto [s4, ds4] = HierAMuS::LobattoShapes::getShape(IntegrationPt.eta, 1);

  auto [s5, ds5] = HierAMuS::LobattoShapes::getShape(IntegrationPt.zeta, 0);
  auto [s6, ds6] = HierAMuS::LobattoShapes::getShape(IntegrationPt.zeta, 1);

  shapes.shapes(0) = s1 * s3 * s5;
  shapes.shapeDeriv(0, 0) = ds1 * s3 * s5;
  shapes.shapeDeriv(1, 0) = s1 * ds3 * s5;
  shapes.shapeDeriv(2, 0) = s1 * s3 * ds5;

  shapes.shapes(1) = s2 * s3 * s5;
  shapes.shapeDeriv(0, 1) = ds2 * s3 * s5;
  shapes.shapeDeriv(1, 1) = s2 * ds3 * s5;
  shapes.shapeDeriv(2, 1) = s2 * s3 * ds5;

  shapes.shapes(2) = s2 * s4 * s5;
  shapes.shapeDeriv(0, 2) = ds2 * s4 * s5;
  shapes.shapeDeriv(1, 2) = s2 * ds4 * s5;
  shapes.shapeDeriv(2, 2) = s2 * s4 * ds5;

  shapes.shapes(3) = s1 * s4 * s5;
  shapes.shapeDeriv(0, 3) = ds1 * s4 * s5;
  shapes.shapeDeriv(1, 3) = s1 * ds4 * s5;
  shapes.shapeDeriv(2, 3) = s1 * s4 * ds5;

  shapes.shapes(4) = s1 * s3 * s6;
  shapes.shapeDeriv(0, 4) = ds1 * s3 * s6;
  shapes.shapeDeriv(1, 4) = s1 * ds3 * s6;
  shapes.shapeDeriv(2, 4) = s1 * s3 * ds6;

  shapes.shapes(5) = s2 * s3 * s6;
  shapes.shapeDeriv(0, 5) = ds2 * s3 * s6;
  shapes.shapeDeriv(1, 5) = s2 * ds3 * s6;
  shapes.shapeDeriv(2, 5) = s2 * s3 * ds6;

  shapes.shapes(6) = s2 * s4 * s6;
  shapes.shapeDeriv(0, 6) = ds2 * s4 * s6;
  shapes.shapeDeriv(1, 6) = s2 * ds4 * s6;
  shapes.shapeDeriv(2, 6) = s2 * s4 * ds6;

  shapes.shapes(7) = s1 * s4 * s6;
  shapes.shapeDeriv(0, 7) = ds1 * s4 * s6;
  shapes.shapeDeriv(1, 7) = s1 * ds4 * s6;
  shapes.shapeDeriv(2, 7) = s1 * s4 * ds6;

  if (order > 1) {
    indexType spos = 8;
    // Edges
    // Edge 1
    IntegrationPoint edgeIP;
    auto geoData = pointers.getGeometryData();
    {
      auto &edge = geoData->getEdge(this->m_edges[0]);
      prec orient = edge.getEdgeOrientation(this->m_verts[0], this->m_verts[1]);
      prec xsiT = IntegrationPt.xi * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge.getH1ShapesInternal(pointers, order, edgeIP);
      for (auto i = 0; i < eShape.shapes.size(); ++i) {
        shapes.shapes(spos) = eShape.shapes(i) * s3 * s5;
        shapes.shapeDeriv(0, spos) = eShape.shapeDeriv(i) * s3 * s5 * orient;
        shapes.shapeDeriv(1, spos) = eShape.shapes(i) * ds3 * s5;
        shapes.shapeDeriv(2, spos) = eShape.shapes(i) * s3 * ds5;
        ++spos;
      }
    }
    // Edge 2
    {
      auto &edge = geoData->getEdge(this->m_edges[1]);
      prec orient = edge.getEdgeOrientation(this->m_verts[1], this->m_verts[2]);
      prec xsiT = IntegrationPt.eta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge.getH1ShapesInternal(pointers, order, edgeIP);
      for (auto i = 0; i < eShape.shapes.size(); ++i) {
        shapes.shapes(spos) = eShape.shapes(i) * s2 * s5;
        shapes.shapeDeriv(0, spos) = eShape.shapes(i) * ds2 * s5;
        shapes.shapeDeriv(1, spos) = eShape.shapeDeriv(i) * s2 * s5 * orient;
        shapes.shapeDeriv(2, spos) = eShape.shapes(i) * s2 * ds5;
        ++spos;
      }
    }
    // Edge 3
    {
      auto &edge = geoData->getEdge(this->m_edges[2]);
      prec orient = edge.getEdgeOrientation(this->m_verts[3], this->m_verts[2]);
      prec xsiT = IntegrationPt.xi * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge.getH1ShapesInternal(pointers, order, edgeIP);
      for (auto i = 0; i < eShape.shapes.size(); ++i) {
        shapes.shapes(spos) = eShape.shapes(i) * s4 * s5;
        shapes.shapeDeriv(0, spos) = eShape.shapeDeriv(i) * s4 * s5 * orient;
        shapes.shapeDeriv(1, spos) = eShape.shapes(i) * ds4 * s5;
        shapes.shapeDeriv(2, spos) = eShape.shapes(i) * s4 * ds5;
        ++spos;
      }
    }
    // Edge 4
    {
      auto &edge = geoData->getEdge(this->m_edges[3]);
      prec orient = edge.getEdgeOrientation(this->m_verts[0], this->m_verts[3]);
      prec xsiT = IntegrationPt.eta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge.getH1ShapesInternal(pointers, order, edgeIP);
      for (auto i = 0; i < eShape.shapes.size(); ++i) {
        shapes.shapes(spos) = eShape.shapes(i) * s1 * s5;
        shapes.shapeDeriv(0, spos) = eShape.shapes(i) * ds1 * s5;
        shapes.shapeDeriv(1, spos) = eShape.shapeDeriv(i) * s1 * s5 * orient;
        shapes.shapeDeriv(2, spos) = eShape.shapes(i) * s1 * ds5;
        ++spos;
      }
    }
    // Edge 5
    {
      auto &edge = geoData->getEdge(this->m_edges[4]);
      prec orient = edge.getEdgeOrientation(this->m_verts[0], this->m_verts[4]);
      prec xsiT = IntegrationPt.zeta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge.getH1ShapesInternal(pointers, order, edgeIP);
      for (auto i = 0; i < eShape.shapes.size(); ++i) {
        shapes.shapes(spos) = eShape.shapes(i) * s1 * s3;
        shapes.shapeDeriv(0, spos) = eShape.shapes(i) * ds1 * s3;
        shapes.shapeDeriv(1, spos) = eShape.shapes(i) * s1 * ds3;
        shapes.shapeDeriv(2, spos) = eShape.shapeDeriv(i) * s1 * s3 * orient;
        ++spos;
      }
    }
    // Edge 6
    {
      auto &edge = geoData->getEdge(this->m_edges[5]);
      prec orient = edge.getEdgeOrientation(this->m_verts[1], this->m_verts[5]);
      prec xsiT = IntegrationPt.zeta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge.getH1ShapesInternal(pointers, order, edgeIP);
      for (auto i = 0; i < eShape.shapes.size(); ++i) {
        shapes.shapes(spos) = eShape.shapes(i) * s2 * s3;
        shapes.shapeDeriv(0, spos) = eShape.shapes(i) * ds2 * s3;
        shapes.shapeDeriv(1, spos) = eShape.shapes(i) * s2 * ds3;
        shapes.shapeDeriv(2, spos) = eShape.shapeDeriv(i) * s2 * s3 * orient;
        ++spos;
      }
    }
    // Edge 7
    {
      auto &edge = geoData->getEdge(this->m_edges[6]);
      prec orient = edge.getEdgeOrientation(this->m_verts[2], this->m_verts[6]);
      prec xsiT = IntegrationPt.zeta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge.getH1ShapesInternal(pointers, order, edgeIP);
      for (auto i = 0; i < eShape.shapes.size(); ++i) {
        shapes.shapes(spos) = eShape.shapes(i) * s2 * s4;
        shapes.shapeDeriv(0, spos) = eShape.shapes(i) * ds2 * s4;
        shapes.shapeDeriv(1, spos) = eShape.shapes(i) * s2 * ds4;
        shapes.shapeDeriv(2, spos) = eShape.shapeDeriv(i) * s2 * s4 * orient;
        ++spos;
      }
    }
    // Edge 8
    {
      auto &edge = geoData->getEdge(this->m_edges[7]);
      prec orient = edge.getEdgeOrientation(this->m_verts[3], this->m_verts[7]);
      prec xsiT = IntegrationPt.zeta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge.getH1ShapesInternal(pointers, order, edgeIP);
      for (auto i = 0; i < eShape.shapes.size(); ++i) {
        shapes.shapes(spos) = eShape.shapes(i) * s1 * s4;
        shapes.shapeDeriv(0, spos) = eShape.shapes(i) * ds1 * s4;
        shapes.shapeDeriv(1, spos) = eShape.shapes(i) * s1 * ds4;
        shapes.shapeDeriv(2, spos) = eShape.shapeDeriv(i) * s1 * s4 * orient;
        ++spos;
      }
    }
    // Edge 9
    {
      auto &edge = geoData->getEdge(this->m_edges[8]);
      prec orient = edge.getEdgeOrientation(this->m_verts[4], this->m_verts[5]);
      prec xsiT = IntegrationPt.xi * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge.getH1ShapesInternal(pointers, order, edgeIP);
      for (auto i = 0; i < eShape.shapes.size(); ++i) {
        shapes.shapes(spos) = eShape.shapes(i) * s3 * s6;
        shapes.shapeDeriv(0, spos) = eShape.shapeDeriv(i) * s3 * s6 * orient;
        shapes.shapeDeriv(1, spos) = eShape.shapes(i) * ds3 * s6;
        shapes.shapeDeriv(2, spos) = eShape.shapes(i) * s3 * ds6;
        ++spos;
      }
    }
    // Edge 10
    {
      auto &edge = geoData->getEdge(this->m_edges[9]);
      prec orient = edge.getEdgeOrientation(this->m_verts[5], this->m_verts[6]);
      prec xsiT = IntegrationPt.eta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge.getH1ShapesInternal(pointers, order, edgeIP);
      for (auto i = 0; i < eShape.shapes.size(); ++i) {
        shapes.shapes(spos) = eShape.shapes(i) * s2 * s6;
        shapes.shapeDeriv(0, spos) = eShape.shapes(i) * ds2 * s6;
        shapes.shapeDeriv(1, spos) = eShape.shapeDeriv(i) * s2 * s6 * orient;
        shapes.shapeDeriv(2, spos) = eShape.shapes(i) * s2 * ds6;
        ++spos;
      }
    }
    // Edge 11
    {
      auto &edge = geoData->getEdge(this->m_edges[10]);
      prec orient = edge.getEdgeOrientation(this->m_verts[7], this->m_verts[6]);
      prec xsiT = IntegrationPt.xi * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge.getH1ShapesInternal(pointers, order, edgeIP);
      for (auto i = 0; i < eShape.shapes.size(); ++i) {
        shapes.shapes(spos) = eShape.shapes(i) * s4 * s6;
        shapes.shapeDeriv(0, spos) = eShape.shapeDeriv(i) * s4 * s6 * orient;
        shapes.shapeDeriv(1, spos) = eShape.shapes(i) * ds4 * s6;
        shapes.shapeDeriv(2, spos) = eShape.shapes(i) * s4 * ds6;
        ++spos;
      }
    }
    // Edge 12
    {
      auto &edge = geoData->getEdge(this->m_edges[11]);
      prec orient = edge.getEdgeOrientation(this->m_verts[4], this->m_verts[7]);
      prec xsiT = IntegrationPt.eta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge.getH1ShapesInternal(pointers, order, edgeIP);
      for (auto i = 0; i < eShape.shapes.size(); ++i) {
        shapes.shapes(spos) = eShape.shapes(i) * s1 * s6;
        shapes.shapeDeriv(0, spos) = eShape.shapes(i) * ds1 * s6;
        shapes.shapeDeriv(1, spos) = eShape.shapeDeriv(i) * s1 * s6 * orient;
        shapes.shapeDeriv(2, spos) = eShape.shapes(i) * s1 * ds6;
        ++spos;
      }
    }

    // Face 1
    {
      auto face = geoData->getFace(this->m_faces[0]);
      auto faceOrtientation =
          face->getOrientation(pointers, this->m_verts[0], this->m_verts[1]);
      IntegrationPoint faceIp;
      faceIp.xi = IntegrationPt.xi;
      faceIp.eta = IntegrationPt.eta;
      auto faceShape =
          face->getH1ShapesInternal(pointers, order, faceIp, faceOrtientation);
      for (auto i = 0; i < faceShape.shapes.size(); ++i) {
        shapes.shapes(spos) = faceShape.shapes(i) * s5;
        shapes.shapeDeriv(0, spos) = faceShape.shapeDeriv(0, i) * s5;
        shapes.shapeDeriv(1, spos) = faceShape.shapeDeriv(1, i) * s5;
        shapes.shapeDeriv(2, spos) = faceShape.shapes(i) * ds5;
        ++spos;
      }
    }
    // Face 2
    {
      auto face = geoData->getFace(this->m_faces[1]);
      auto faceOrtientation =
          face->getOrientation(pointers, this->m_verts[0], this->m_verts[4]);
      IntegrationPoint faceIp;
      faceIp.xi = IntegrationPt.zeta;
      faceIp.eta = IntegrationPt.xi;
      auto faceShape =
          face->getH1ShapesInternal(pointers, order, faceIp, faceOrtientation);
      for (auto i = 0; i < faceShape.shapes.size(); ++i) {
        shapes.shapes(spos) = faceShape.shapes(i) * s3;
        shapes.shapeDeriv(0, spos) = faceShape.shapeDeriv(1, i) * s3;
        shapes.shapeDeriv(1, spos) = faceShape.shapes(i) * ds3;
        shapes.shapeDeriv(2, spos) = faceShape.shapeDeriv(0, i) * s3;
        ++spos;
      }
    }
    // Face 3
    {
      auto face = geoData->getFace(this->m_faces[2]);
      auto faceOrtientation =
          face->getOrientation(pointers, this->m_verts[1], this->m_verts[2]);
      IntegrationPoint faceIp;
      faceIp.xi = IntegrationPt.eta;
      faceIp.eta = IntegrationPt.zeta;
      auto faceShape =
          face->getH1ShapesInternal(pointers, order, faceIp, faceOrtientation);
      for (auto i = 0; i < faceShape.shapes.size(); ++i) {
        shapes.shapes(spos) = faceShape.shapes(i) * s2;
        shapes.shapeDeriv(0, spos) = faceShape.shapes(i) * ds2;
        shapes.shapeDeriv(1, spos) = faceShape.shapeDeriv(0, i) * s2;
        shapes.shapeDeriv(2, spos) = faceShape.shapeDeriv(1, i) * s2;
        ++spos;
      }
    }
    // Face 4
    {
      auto face = geoData->getFace(this->m_faces[3]);
      auto faceOrtientation =
          face->getOrientation(pointers, this->m_verts[3], this->m_verts[7]);
      IntegrationPoint faceIp;
      faceIp.xi = IntegrationPt.zeta;
      faceIp.eta = IntegrationPt.xi;
      auto faceShape =
          face->getH1ShapesInternal(pointers, order, faceIp, faceOrtientation);
      for (auto i = 0; i < faceShape.shapes.size(); ++i) {
        shapes.shapes(spos) = faceShape.shapes(i) * s4;
        shapes.shapeDeriv(0, spos) = faceShape.shapeDeriv(1, i) * s4;
        shapes.shapeDeriv(1, spos) = faceShape.shapes(i) * ds4;
        shapes.shapeDeriv(2, spos) = faceShape.shapeDeriv(0, i) * s4;
        ++spos;
      }
    }
    // Face 5
    {
      auto face = geoData->getFace(this->m_faces[4]);
      auto faceOrtientation =
          face->getOrientation(pointers, this->m_verts[0], this->m_verts[3]);
      IntegrationPoint faceIp;
      faceIp.xi = IntegrationPt.eta;
      faceIp.eta = IntegrationPt.zeta;
      auto faceShape =
          face->getH1ShapesInternal(pointers, order, faceIp, faceOrtientation);
      for (auto i = 0; i < faceShape.shapes.size(); ++i) {
        shapes.shapes(spos) = faceShape.shapes(i) * s1;
        shapes.shapeDeriv(0, spos) = faceShape.shapes(i) * ds1;
        shapes.shapeDeriv(1, spos) = faceShape.shapeDeriv(0, i) * s1;
        shapes.shapeDeriv(2, spos) = faceShape.shapeDeriv(1, i) * s1;
        ++spos;
      }
    }
    // Face 6
    {
      auto face = geoData->getFace(this->m_faces[5]);
      auto faceOrtientation =
          face->getOrientation(pointers, this->m_verts[4], this->m_verts[5]);
      IntegrationPoint faceIp;
      faceIp.xi = IntegrationPt.xi;
      faceIp.eta = IntegrationPt.eta;
      auto faceShape =
          face->getH1ShapesInternal(pointers, order, faceIp, faceOrtientation);
      for (auto i = 0; i < faceShape.shapes.size(); ++i) {
        shapes.shapes(spos) = faceShape.shapes(i) * s6;
        shapes.shapeDeriv(0, spos) = faceShape.shapeDeriv(0, i) * s6;
        shapes.shapeDeriv(1, spos) = faceShape.shapeDeriv(1, i) * s6;
        shapes.shapeDeriv(2, spos) = faceShape.shapes(i) * ds6;
        ++spos;
      }
    }
    auto internalShapes =
        this->getH1ShapesInternal(pointers, order, IntegrationPt);
    for (auto i = 0; i < internalShapes.shapes.size(); ++i) {
      shapes.shapes(spos) = internalShapes.shapes(i);
      shapes.shapeDeriv(0, spos) = internalShapes.shapeDeriv(0, i);
      shapes.shapeDeriv(1, spos) = internalShapes.shapeDeriv(1, i);
      shapes.shapeDeriv(2, spos) = internalShapes.shapeDeriv(2, i);
      ++spos;
    }
  }

  return shapes;
};

auto LinearBrick::getH1ShapesInternal(PointerCollection &pointers,
                                      indexType order,
                                      IntegrationPoint &IntegrationPt)
    -> H1Shapes {
  H1Shapes bubbleShapes;
  indexType numNodes = order - 1;
  indexType numShapes = numNodes * numNodes * numNodes;
  bubbleShapes.shapes.resize(numShapes);
  bubbleShapes.shapeDeriv.resize(3, numShapes);
  indexType pos = 0;
  for (auto i = 2; i < numNodes + 2; ++i) {
    auto xiS = LobattoShapes::getShape(IntegrationPt.xi, i);
    for (auto j = 2; j < numNodes + 2; ++j) {
      auto etaS = LobattoShapes::getShape(IntegrationPt.eta, j);
      for (auto k = 2; k < numNodes + 2; ++k) {
        auto zetaS = LobattoShapes::getShape(IntegrationPt.zeta, k);
        bubbleShapes.shapes(pos) =
            xiS.shapeValue * etaS.shapeValue * zetaS.shapeValue;
        bubbleShapes.shapeDeriv(0, pos) =
            xiS.shapeDerivative * etaS.shapeValue * zetaS.shapeValue;
        bubbleShapes.shapeDeriv(1, pos) =
            xiS.shapeValue * etaS.shapeDerivative * zetaS.shapeValue;
        bubbleShapes.shapeDeriv(2, pos) =
            xiS.shapeValue * etaS.shapeValue * zetaS.shapeDerivative;
        ++pos;
      }
    }
  }
  return bubbleShapes;
}

auto LinearBrick::getH1Nodes(PointerCollection &pointers, indexType meshID,
                             indexType order) -> std::vector<GenericNodes *> {
  return {};
}

auto LinearBrick::getH1NodesInternal(PointerCollection &pointers,
                                     indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  if (order > 1) {
    indexType totnodes = (order - 1) * (order - 1) * (order - 1);
    auto tempnodes = this->getNodesOfSet(pointers, meshID);
    nodes.insert(nodes.end(), tempnodes.begin(), tempnodes.end());
    if (totnodes != tempnodes.size()) {
      std::stringstream ss;
      ss << "Error in LinearQuadrilateral::getH1NodesInternal: "
         << "Expected " << totnodes << " nodes, got " << tempnodes.size()
         << " instead.";
      throw std::runtime_error(ss.str());
    }
  }
  return nodes;
}

void LinearBrick::getFaces(PointerCollection &pointers,
                           std::vector<Base *> &facesOut) {
  facesOut.clear();
  for (auto i = 0; i < 12; ++i) {
    facesOut.push_back(&pointers.getGeometryData()->getEdge(this->m_edges[i]));
  }
}

void LinearBrick::geometryToParaview(PointerCollection &pointers,
                                     vtkPlotInterface &paraviewAdapter,
                                     indexType mainMesh, indexType subMesh) {
  indexType numPoints = 8;
  std::vector<indexType> points(numPoints);
  points.clear();
  for (auto i : this->m_verts) {
    auto &vert = pointers.getGeometryData()->getVertex(i);
    vert.geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
    points.push_back(i);
  }

  paraviewAdapter.addCell(mainMesh, subMesh, this->id, 1, points, numPoints,
                          VTK_HEXAHEDRON);
}

void LinearBrick::computeWeightsParaview(PointerCollection &pointers,
                                         vtkPlotInterface &paraviewAdapter,
                                         indexType mainMesh,
                                         indexType subMesh) {
  auto GP = this->getIntegrationPoints(pointers, -1);
  GP.setOrder(2);

  for (auto i : GP) {
    auto jaco = this->getJacobian(pointers, i);
    auto shapes = this->getH1Shapes(pointers, 1, i);
    prec dA = jaco.determinant() * i.weight;

    for (auto i = 0; i < 8; ++i) {
      auto &vert = pointers.getGeometryData()->getVertex(this->m_verts[i]);
      std::vector<prec> val;
      val.push_back(shapes.shapes(i) * dA);
      paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, val, vert.getId(),
                                           1, paraviewNames::weightName());
    }
  }
}

void LinearBrick::H1SolutionToParaview(PointerCollection &pointers,
                                       vtkPlotInterface &paraviewAdapter,
                                       indexType mainMesh, indexType subMesh,
                                       indexType meshId, indexType order,
                                       std::string &name) {
  std::vector<DegreeOfFreedom *> Dofs;
  this->getH1Dofs(pointers, Dofs, meshId, order);
  auto solution = pointers.getSolutionState()->getSolution(Dofs);
  for (auto i = 0; i < 8; ++i) {
    auto &V = pointers.getGeometryData()->getVertex(this->m_verts[i]);
    std::vector<prec> sol(3);
    for (auto j = 0; j < 3; ++j) {
      sol[j] = solution(3 * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, V.getId(), sol, 3, name);
  }
}

void LinearBrick::H1DataToParaview(PointerCollection &pointers,
                                   vtkPlotInterface &paraviewAdapter,
                                   indexType mainMesh, indexType subMesh,
                                   Types::VectorX<prec> &Data,
                                   indexType numberComponents, indexType order,
                                   std::string &name) {
  for (auto i = 0; i < 8; ++i) {
    auto &V = pointers.getGeometryData()->getVertex(this->m_verts[i]);
    std::vector<prec> sol(numberComponents);
    for (auto j = 0; j < numberComponents; ++j) {
      sol[j] = Data(numberComponents * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, V.getId(), sol,
                                 numberComponents, name);
  }
}

void LinearBrick::projectDataToParaviewVertices(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {

  auto shapes = this->getH1Shapes(pointers, 1, IntegrationPt);
  std::vector<prec> vals(numberComponents);

  auto jaco = this->getJacobian(pointers, IntegrationPt);
  auto dA = jaco.determinant() * IntegrationPt.weight;
  for (auto i = 0; i < 8; ++i) {
    auto &V = pointers.getGeometryData()->getVertex(this->m_verts[i]);
    for (auto j = 0; j < numberComponents; ++j) {
      vals[j] = data(j) * shapes.shapes(i) * dA;
    }
    paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, vals, V.getId(),
                                         numberComponents, name);
  }
}

void LinearBrick::checkUpdateElement(GeometryData &geoData) {

  static const std::vector<indexType> lv1({0, 0, 1, 2, 3, 4});
  static const std::vector<indexType> lv2({1, 1, 2, 3, 0, 5});
  static const std::vector<indexType> lv3({2, 5, 6, 7, 4, 6});
  static const std::vector<indexType> lv4({3, 4, 5, 6, 7, 7});

  for (auto i = 0; i < 6; ++i) {
    if (this->m_faces[i] == -1) {
      auto &V1 = geoData.getVertex(this->m_verts[lv1[i]]);
      auto &V2 = geoData.getVertex(this->m_verts[lv2[i]]);
      auto &V3 = geoData.getVertex(this->m_verts[lv3[i]]);
      auto &V4 = geoData.getVertex(this->m_verts[lv4[i]]);

      auto faceNums = V1.getConnectedFaces();
      bool faceExists = false;
      bool search = true;
      if (faceNums.empty()) {
        search = false;
      }
      indexType pos = 0;
      while (search) {
        auto face = geoData.getFace(faceNums[pos]);
        if (face->hasVertices(V1.getId(), V2.getId(), V3.getId())) {
          faceExists = true;
          search = false;
          this->m_faces[i] = face->getId();
        } else {
          ++pos;
          if (pos >= faceNums.size()) {
            search = false;
          }
        }
      }
      if (!faceExists) {
        auto facenum = geoData.requestNewGeometryObject(
            GeometryTypes::LinearQuadrilateral);
        auto face = geoData.getFace(facenum);
        std::vector<indexType> verts(
            {V1.getId(), V2.getId(), V3.getId(), V4.getId()});
        face->setVerts(geoData, verts);
        this->m_faces[i] = face->getId();
        face->checkUpdateElement(geoData);
      }
    }
  }

  // Edges
  static const std::vector<indexType> le1({0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7});
  static const std::vector<indexType> le2({1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 7, 4});
  for (auto i = 0; i < 12; ++i) {
    if (this->m_edges[i] == -1) {
      auto &V1 = geoData.getVertex(this->m_verts[le1[i]]);
      auto &V2 = geoData.getVertex(this->m_verts[le2[i]]);

      bool search = true;
      bool edgeExists = false;
      auto edgeNums = V1.getConnectedEdges();
      if (edgeNums.empty()) {
        search = false;
      }
      indexType pos = 0;
      while (search) {
        auto &edge = geoData.getEdge(edgeNums[pos]);
        if (edge.hasVertices(V1.getId(), V2.getId())) {
          edgeExists = true;
          search = false;
          this->m_edges[i] = edge.getId();
        } else {
          pos++;
          if (pos >= edgeNums.size()) {
            search = false;
          }
        }
      }
      if (!edgeExists) {
        std::cout << "something went wrong" << std::endl;
        auto edgeNum =
            geoData.requestNewGeometryObject(GeometryTypes::LinearEdge);
        auto &edge = geoData.getEdge(edgeNum);
        std::vector<indexType> verts = {V1.getId(), V2.getId()};
        edge.setVerts(geoData, verts);
        this->m_edges[i] = edgeNum;
      }
    }
  }
}

const GeometryTypes LinearBrick::type = GeometryTypes::LinearBrick;
} // namespace HierAMuS::Geometry
