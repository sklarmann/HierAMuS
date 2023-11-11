// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "geometry/GeometryBaseData.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

#include <geometry/GeometryData.h>

#include <geometry/Edges/EdgesData.h>
#include <geometry/Faces/FacesData.h>
#include <geometry/VertexData.h>

#include "geometry/Volumes/LinearBrickData.h"
#include "geometry/Volumes/LinearBrickRuntime.h"

#include "geometry/Edges/EdgeH1ShapesInterface.h"
#include "geometry/Edges/EdgesRuntime.h"
#include "geometry/Faces/FacesH1Interface.h"
#include "geometry/Faces/FacesRuntime.h"

#include <types/MatrixTypes.h>

#include <vector>

#include <iomanip>

#include "plot/vtkplotClass.h"

#include "shapefunctions/LobattoShapes.h"
#include "shapefunctions/IntegrationsPoints/helperClasses/IntegrationPointsManagement.h"

#include <vtkCellType.h>
namespace HierAMuS::Geometry {

LinearBrickRuntime::LinearBrickRuntime(GeometryData &geoData,
                                       LinearBrickData &data_element)
    : VolumesRuntimeDataInterface(geoData, data_element),
      m_LinearBrick_Data(data_element) {}

LinearBrickRuntime::~LinearBrickRuntime() = default;

auto LinearBrickRuntime::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::LinearBrickRuntime::type;
}

void LinearBrickRuntime::print(spdlog::logger &Log) {
  m_LinearBrick_Data.print(Log);
}

void LinearBrickRuntime::getEdgeNumbers(std::vector<indexType> &edgesOut) {
  m_LinearBrick_Data.getEdgeNumbers(edgesOut);
}

void LinearBrickRuntime::getFaceNumbers(std::vector<indexType> &facesOut) {
  m_LinearBrick_Data.getFaceNumbers(facesOut);
}

auto LinearBrickRuntime::getIntegrationPoints(indexType elementId)
    -> IntegrationPoints {
  auto points = IntegrationPointsManagement::getIntegrationsPoints(elementId);
  points.setType(IntegrationType::Gauss3D);
  return points;
}

void LinearBrickRuntime::setH1Shapes(indexType meshId, indexType order,
                                     NodeTypes type) {
  for (auto i = 0; i < m_numberOfFaces; i++) {
    m_Faces[i]->getH1Face()->setH1Shapes(meshId, order, type);
  }
  if (order > 1) {
    this->setH1ShapesInternal(meshId, order, type);
  }
}

void LinearBrickRuntime::setH1ShapesInternal(indexType meshId, indexType order,
                                             NodeTypes type) {
  if (order > 1) {
    indexType num = order - 1;
    num *= num * num;
    this->setNodeSet(meshId, num, type);
  }
}

void LinearBrickRuntime::getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs,
                                   indexType meshID, indexType order) {

  indexType pp = order + 1;
  pp *= pp * pp * 3;
  Dofs.reserve(pp);

  for (auto i = 0; i < m_numberOfVerts; i++) {
    auto nodeList = m_Vertices[i]->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }

  if (order > 1) {
    for (auto i = 0; i < 12; i++) {
      m_Edges[i]->getH1Edge()->getH1DofsInternal(Dofs, meshID, order);
    }

    for (auto i = 0; i < m_numberOfFaces; i++) {
      m_Faces[i]->getH1Face()->getH1DofsInternal(Dofs, meshID, order);
    }
    this->getH1DofsInternal(Dofs, meshID, order);
  }
}

void LinearBrickRuntime::getH1DofsInternal(std::vector<DegreeOfFreedom *> &Dofs,
                                           indexType meshID, indexType order) {
  if (order > 1) {
    auto nodeList = this->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }
}

void LinearBrickRuntime::getH1Shapes(indexType order,
                                     Types::VectorX<prec> &shape,
                                     Types::Matrix3X<prec> &shapeDerivative,
                                     prec xsi, prec eta, prec zeta) {
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
    this->getH1ShapesInternal(order, tempshape, tempdshape, xsi, eta,
                              zeta);
  }
}

void LinearBrickRuntime::getH1ShapesInternal(
    indexType order, Types::VectorX<prec> &shape,
    Types::Matrix3X<prec> &shapeDerivative, prec xsi, prec eta, prec zeta) {
  if (order > 1) {
    throw std::runtime_error(
        "Higher order shapes not implemented for Brick element!");
  }
}

auto LinearBrickRuntime::getH1Shapes(indexType order,
                                     IntegrationPoint &IntegrationPt)
    -> H1Shapes {

  indexType numshapes = (order + 1);
  numshapes *= numshapes * numshapes;

  H1Shapes shapes(numshapes, 3);

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
    {
      auto edge = m_Edges[0];
      prec orient = edge->getEdgeOrientation(this->m_Vertices[0]->getId(),
                                             this->m_Vertices[1]->getId());
      prec xsiT = IntegrationPt.xi * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge->getH1Edge()->getH1ShapesInternal(order, edgeIP);
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
      auto edge = m_Edges[1];
      prec orient = edge->getEdgeOrientation(this->m_Vertices[1]->getId(),
                                             this->m_Vertices[2]->getId());
      prec xsiT = IntegrationPt.eta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge->getH1Edge()->getH1ShapesInternal(order, edgeIP);
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
      auto edge = m_Edges[2];
      prec orient = edge->getEdgeOrientation(this->m_Vertices[3]->getId(),
                                             this->m_Vertices[2]->getId());
      prec xsiT = IntegrationPt.xi * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge->getH1Edge()->getH1ShapesInternal(order, edgeIP);
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
      auto edge = m_Edges[3];
      prec orient = edge->getEdgeOrientation(this->m_Vertices[0]->getId(),
                                             this->m_Vertices[3]->getId());
      prec xsiT = IntegrationPt.eta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge->getH1Edge()->getH1ShapesInternal(order, edgeIP);
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
      auto &edge = m_Edges[4];
      prec orient = edge->getEdgeOrientation(this->m_Vertices[0]->getId(),
                                             this->m_Vertices[4]->getId());
      prec xsiT = IntegrationPt.zeta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge->getH1Edge()->getH1ShapesInternal(order, edgeIP);
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
      auto edge = m_Edges[5];
      prec orient = edge->getEdgeOrientation(this->m_Vertices[1]->getId(),
                                             this->m_Vertices[5]->getId());
      prec xsiT = IntegrationPt.zeta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge->getH1Edge()->getH1ShapesInternal(order, edgeIP);
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
      auto edge = m_Edges[6];
      prec orient = edge->getEdgeOrientation(this->m_Vertices[2]->getId(),
                                             this->m_Vertices[6]->getId());
      prec xsiT = IntegrationPt.zeta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge->getH1Edge()->getH1ShapesInternal(order, edgeIP);
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
      auto edge = m_Edges[7];
      prec orient = edge->getEdgeOrientation(this->m_Vertices[3]->getId(),
                                             this->m_Vertices[7]->getId());
      prec xsiT = IntegrationPt.zeta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge->getH1Edge()->getH1ShapesInternal(order, edgeIP);
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
      auto edge = m_Edges[8];
      prec orient = edge->getEdgeOrientation(this->m_Vertices[4]->getId(),
                                             this->m_Vertices[5]->getId());
      prec xsiT = IntegrationPt.xi * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge->getH1Edge()->getH1ShapesInternal(order, edgeIP);
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
      auto edge = m_Edges[9];
      prec orient = edge->getEdgeOrientation(this->m_Vertices[5]->getId(),
                                             this->m_Vertices[6]->getId());
      prec xsiT = IntegrationPt.eta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge->getH1Edge()->getH1ShapesInternal(order, edgeIP);
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
      auto edge = m_Edges[10];
      prec orient = edge->getEdgeOrientation(this->m_Vertices[7]->getId(),
                                             this->m_Vertices[6]->getId());
      prec xsiT = IntegrationPt.xi * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge->getH1Edge()->getH1ShapesInternal(order, edgeIP);
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
      auto edge = m_Edges[11];
      prec orient = edge->getEdgeOrientation(this->m_Vertices[4]->getId(),
                                             this->m_Vertices[7]->getId());
      prec xsiT = IntegrationPt.eta * orient;
      edgeIP.xi = xsiT;
      auto eShape = edge->getH1Edge()->getH1ShapesInternal(order, edgeIP);
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
      auto face = m_Faces[0];
      auto faceOrtientation = face->getOrientation(
          this->m_Vertices[0]->getId(), this->m_Vertices[1]->getId());
      IntegrationPoint faceIp;
      faceIp.xi = IntegrationPt.xi;
      faceIp.eta = IntegrationPt.eta;
      auto faceShape = face->getH1Face()->getH1ShapesInternal(order, faceIp,
                                                              faceOrtientation);
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
      auto face = m_Faces[1];
      auto faceOrtientation = face->getOrientation(
          this->m_Vertices[0]->getId(), this->m_Vertices[4]->getId());
      IntegrationPoint faceIp;
      faceIp.xi = IntegrationPt.zeta;
      faceIp.eta = IntegrationPt.xi;
      auto faceShape = face->getH1Face()->getH1ShapesInternal(order, faceIp,
                                                              faceOrtientation);
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
      auto face = m_Faces[2];
      auto faceOrtientation = face->getOrientation(
          this->m_Vertices[1]->getId(), this->m_Vertices[2]->getId());
      IntegrationPoint faceIp;
      faceIp.xi = IntegrationPt.eta;
      faceIp.eta = IntegrationPt.zeta;
      auto faceShape = face->getH1Face()->getH1ShapesInternal(order, faceIp,
                                                              faceOrtientation);
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
      auto face = m_Faces[3];
      auto faceOrtientation = face->getOrientation(
          this->m_Vertices[3]->getId(), this->m_Vertices[7]->getId());
      IntegrationPoint faceIp;
      faceIp.xi = IntegrationPt.zeta;
      faceIp.eta = IntegrationPt.xi;
      auto faceShape = face->getH1Face()->getH1ShapesInternal(order, faceIp,
                                                              faceOrtientation);
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
      auto face = m_Faces[4];
      auto faceOrtientation = face->getOrientation(
          this->m_Vertices[0]->getId(), this->m_Vertices[3]->getId());
      IntegrationPoint faceIp;
      faceIp.xi = IntegrationPt.eta;
      faceIp.eta = IntegrationPt.zeta;
      auto faceShape = face->getH1Face()->getH1ShapesInternal(order, faceIp,
                                                              faceOrtientation);
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
      auto face = m_Faces[5];
      auto faceOrtientation = face->getOrientation(
          this->m_Vertices[4]->getId(), this->m_Vertices[5]->getId());
      IntegrationPoint faceIp;
      faceIp.xi = IntegrationPt.xi;
      faceIp.eta = IntegrationPt.eta;
      auto faceShape = face->getH1Face()->getH1ShapesInternal(order, faceIp,
                                                              faceOrtientation);
      for (auto i = 0; i < faceShape.shapes.size(); ++i) {
        shapes.shapes(spos) = faceShape.shapes(i) * s6;
        shapes.shapeDeriv(0, spos) = faceShape.shapeDeriv(0, i) * s6;
        shapes.shapeDeriv(1, spos) = faceShape.shapeDeriv(1, i) * s6;
        shapes.shapeDeriv(2, spos) = faceShape.shapes(i) * ds6;
        ++spos;
      }
    }
    auto internalShapes =
        this->getH1ShapesInternal(order, IntegrationPt);
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

auto LinearBrickRuntime::getH1ShapesInternal(indexType order,
                                             IntegrationPoint &IntegrationPt)
    -> H1Shapes {
  indexType numNodes = order - 1;
  indexType numShapes = numNodes * numNodes * numNodes;
  H1Shapes bubbleShapes(numShapes, 3);
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

auto LinearBrickRuntime::getH1Nodes(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  return {};
}

auto LinearBrickRuntime::getH1NodesInternal(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  if (order > 1) {
    indexType totnodes = (order - 1) * (order - 1) * (order - 1);
    auto tempnodes = this->getNodesOfSet(meshID);
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

void LinearBrickRuntime::geometryToParaview(vtkPlotInterface &paraviewAdapter,
                                            indexType mainMesh,
                                            indexType subMesh) {
  indexType numPoints = 8;
  std::vector<indexType> points(numPoints);
  points.clear();
  for (auto &i : this->m_Vertices) {
    i->geometryToParaview(paraviewAdapter, mainMesh, subMesh);
    points.push_back(i->getId());
  }

  paraviewAdapter.addCell(mainMesh, subMesh, m_LinearBrick_Data.getId(), 1,
                          points, numPoints, VTK_HEXAHEDRON);
}

void LinearBrickRuntime::computeWeightsParaview(
    vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh) {
  auto GP = this->getIntegrationPoints(-1);
  GP.setOrder(2);

  for (auto i : GP) {
    auto jaco = this->getJacobian(i);
    auto shapes = this->getH1Shapes(1, i);
    prec dA = jaco.determinant() * i.weight;

    for (auto i = 0; i < 8; ++i) {
      std::vector<prec> val;
      val.push_back(shapes.shapes(i) * dA);
      paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, val,
                                           m_Vertices[i]->getId(), 1,
                                           paraviewNames::weightName());
    }
  }
}

void LinearBrickRuntime::H1SolutionToParaview(
    vtkPlotInterface &paraviewAdapter, indexType mainMesh, indexType subMesh,
    indexType order, Types::VectorX<prec> &solution, std::string &name) {
  for (auto i = 0; i < 8; ++i) {
    std::vector<prec> sol(3);
    for (auto j = 0; j < 3; ++j) {
      sol[j] = solution(3 * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, m_Vertices[i]->getId(), sol,
                                 3, name);
  }
}

void LinearBrickRuntime::H1DataToParaview(vtkPlotInterface &paraviewAdapter,
                                          indexType mainMesh, indexType subMesh,
                                          Types::VectorX<prec> &Data,
                                          indexType numberComponents,
                                          indexType order, std::string &name) {
  for (auto i = 0; i < 8; ++i) {
    std::vector<prec> sol(numberComponents);
    for (auto j = 0; j < numberComponents; ++j) {
      sol[j] = Data(numberComponents * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, m_Vertices[i]->getId(), sol,
                                 numberComponents, name);
  }
}

void LinearBrickRuntime::projectDataToParaviewVertices(
    vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {

  auto shapes = this->getH1Shapes(1, IntegrationPt);
  std::vector<prec> vals(numberComponents);

  auto jaco = this->getJacobian(IntegrationPt);
  auto dA = jaco.determinant() * IntegrationPt.weight;
  for (auto i = 0; i < 8; ++i) {
    for (auto j = 0; j < numberComponents; ++j) {
      vals[j] = data(j) * shapes.shapes(i) * dA;
    }
    paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, vals,
                                         m_Vertices[i]->getId(),
                                         numberComponents, name);
  }
}

const GeometryTypes LinearBrickRuntime::type = GeometryTypes::LinearBrick;
} // namespace HierAMuS::Geometry
