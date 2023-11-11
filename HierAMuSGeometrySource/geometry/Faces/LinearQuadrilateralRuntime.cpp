// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "geometry/GeometryData.h"
#include "MatrixTypes.h"
#include "geometry/Faces/FacesData.h"
#include "geometry/GeometryBaseData.h"
#include "geometry/GeometryTypes.h"
#include "plot/vtkplotClass.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "shapefunctions/LegendreShapes.h"
#include "shapefunctions/LobattoShapes.h"
#include <exception>


#include "GenericNodes.h"

#include <geometry/Faces/LinearQuadrilateralRuntime.h>
#include <sstream>

#include <geometry/Edges/EdgesData.h>
#include "geometry/Edges/EdgesRuntime.h"
#include "geometry/Edges/EdgeH1ShapesInterface.h"
#include <geometry/VertexData.h>
#include "geometry/Faces/LinearQuadrilateralData.h"


#include "LoadList.h"

#include <stdexcept>
#include <vector>
#include <vtkCellType.h>

#include "HelperFunctions.h"

#include "shapefunctions/IntegrationsPoints/helperClasses/IntegrationPointsManagement.h"

namespace HierAMuS::Geometry {

LinearQuadrilateralRuntime::LinearQuadrilateralRuntime(
    GeometryData &geoData, LinearQuadrilateralData &base_element)
    : FacesRuntimeDataInterface(geoData, base_element) {};

LinearQuadrilateralRuntime::~LinearQuadrilateralRuntime() = default;

auto LinearQuadrilateralRuntime::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::LinearQuadrilateralRuntime::type;
}



void LinearQuadrilateralRuntime::print(spdlog::logger &Log) {
  m_Face_Data_Element.print(Log);
}



auto LinearQuadrilateralRuntime::getTangent_G1(
    IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  // auto shapes = this->getH1Shapes(pointers, 1, integrationPoint);
  Types::Vector3<prec> ret;

  auto [x1, dx1] = LobattoShapes::getShape(integrationPoint.xi, 0);
  auto [y1, dy1] = LobattoShapes::getShape(integrationPoint.eta, 0);
  auto [x2, dx2] = LobattoShapes::getShape(integrationPoint.xi, 1);
  auto [y2, dy2] = LobattoShapes::getShape(integrationPoint.eta, 1);

  ret = m_Vertices[0]->getCoordinates() * (dx1 * y1);
  ret += m_Vertices[1]->getCoordinates() * (dx2 * y1);
  ret += m_Vertices[2]->getCoordinates() * (dx2 * y2);
  ret += m_Vertices[3]->getCoordinates() * (dx1 * y2);

  return ret;
}

auto LinearQuadrilateralRuntime::getTangent_G2(
    IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  // auto shapes = this->getH1Shapes(pointers, 1, integrationPoint);
  Types::Vector3<prec> ret;

  auto [x1, dx1] = LobattoShapes::getShape(integrationPoint.xi, 0);
  auto [y1, dy1] = LobattoShapes::getShape(integrationPoint.eta, 0);
  auto [x2, dx2] = LobattoShapes::getShape(integrationPoint.xi, 1);
  auto [y2, dy2] = LobattoShapes::getShape(integrationPoint.eta, 1);

  ret = m_Vertices[0]->getCoordinates() * (x1 * dy1);
  ret += m_Vertices[1]->getCoordinates() * (x2 * dy1);
  ret += m_Vertices[2]->getCoordinates() * (x2 * dy2);
  ret += m_Vertices[3]->getCoordinates() * (x1 * dy2);
  return ret;
}

auto LinearQuadrilateralRuntime::getFaceNormal()
    -> Types::Vector3<prec> {
  auto &V1 = m_Vertices[0];
  auto &V2 = m_Vertices[1];
  auto &V3 = m_Vertices[3];

  Types::Vector3<prec> dx = V2->getCoordinates() - V1->getCoordinates();
  Types::Vector3<prec> dy = V3->getCoordinates() - V1->getCoordinates();

  Types::Vector3<prec> n = dx.cross(dy).normalized();
  return n;
}

auto LinearQuadrilateralRuntime::getOrientation(indexType vertex1,
                                                indexType vertex2)
    -> faceorientation {
  return m_Face_Data_Element.getOrientation(vertex1, vertex2);
}

void LinearQuadrilateralRuntime::modifyIntegrationpoint(
    IntegrationPoint &IP, prec &shapeFactor, faceorientation orientation) {
  shapeFactor = prec(1);
  if (orientation == faceorientation::p_2) {
    prec teta = IP.eta;
    IP.eta = -IP.xi;
    IP.xi = teta;
  } else if (orientation == faceorientation::p_3) {
    IP.xi = -IP.xi;
    IP.eta = -IP.eta;
  } else if (orientation == faceorientation::p_4) {
    IP.eta = IP.eta - shapeFactor;
  } else if (orientation == faceorientation::n_1) {
    IP.xi = IP.xi - shapeFactor;
  } else if (orientation == faceorientation::n_2) {
    IP.eta = IP.eta - shapeFactor;
  } else if (orientation == faceorientation::n_3) {
    IP.xi = IP.xi + shapeFactor;
  } else if (orientation == faceorientation::n_4) {
    IP.eta = IP.eta + shapeFactor;
  }
}

auto LinearQuadrilateralRuntime::getIntegrationPoints(
    indexType elementId) -> IntegrationPoints {
  IntegrationPoints temp =
      IntegrationPointsManagement::getIntegrationsPoints(elementId);
  temp.setType(IntegrationType::Gauss2D);
  return temp;
}


// H1Shapes

void LinearQuadrilateralRuntime::setH1Shapes(indexType meshId, indexType order,
                                             NodeTypes type) {
  for (auto i = 0; i < 4; ++i) {
    m_Edges[i]->getH1Edge()->setH1Shapes(meshId, order, type);
  }
  this->setH1ShapesInternal(meshId, order, type);
}

void LinearQuadrilateralRuntime::setH1ShapesInternal(
    indexType meshId, indexType order,
    NodeTypes type) {
  if (order > 1) {
    indexType numNodes = order - 1;
    numNodes *= numNodes;

    this->setNodeSet(meshId, numNodes, type);
  }
}

auto LinearQuadrilateralRuntime::getH1Dofs(indexType meshID, indexType order)
    -> std::vector<DegreeOfFreedom *> {
  // std::vector<DegreeOfFreedom *> Dofs;
  // this->getH1Dofs(pointers, Dofs, meshID, order);
  // return Dofs;
  auto nodeList = this->getH1NodesList(meshID, order);
  return nodeList.getDegreesOfFreedom();
}

void LinearQuadrilateralRuntime::getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs,
                                           indexType meshID, indexType order) {
  indexType pp = order + 1;
  pp *= pp * 3;
  Dofs.reserve(pp);
  for (auto i = 0; i < 4; ++i) {
    auto nodeList = m_Vertices[i]->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }
  if (order > 1) {
    for (auto i = 0; i < 4; ++i) {
      m_Edges[i]->getH1Edge()->getH1DofsInternal(Dofs, meshID, order);
    }
    this->getH1DofsInternal(Dofs, meshID, order);
  }
}

void LinearQuadrilateralRuntime::getH1DofsInternal(
    std::vector<DegreeOfFreedom *> &Dofs,
    indexType meshID, indexType order) {
  if (order > 1) {
    auto nodeList = this->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }
}

auto LinearQuadrilateralRuntime::getH1Nodes(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  for (auto i = 0; i < 4; ++i) {
    auto tempNodes = m_Vertices[i]->getNodesOfSet(meshID);
    nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
  }
  if (order > 1) {
    for (auto i = 0; i < 4; ++i) {
      auto tempNodes =
          m_Edges[i]->getH1Edge()->getH1NodesInternal(meshID, order);
      nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
    }
    auto tempNodes = this->getH1NodesInternal(meshID, order);
    nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
  }
  return nodes;
}

auto LinearQuadrilateralRuntime::getH1NodesInternal(indexType meshID,
                                                    indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  if (order > 1) {
    indexType totnodes = (order - 1) * (order - 1);
    auto nodeList = this->getNodeSetNodeListMeshId(meshID);
    nodes = nodeList.getNodes();

    if (totnodes != nodes.size()) {
      std::stringstream ss;
      ss << "Error in LinearQuadrilateral::getH1NodesInternal: "
         << "Expected " << totnodes << " nodes, got " << nodes.size()
         << " instead.";
      throw std::runtime_error(ss.str());
    }
  }
  return nodes;
}

auto LinearQuadrilateralRuntime::getH1NodesList(indexType meshID,
                                                indexType order)
    -> MeshIdNodeList {

  MeshIdNodeList nodeList(meshID);
  nodeList.reserve(9);

  for (auto& i : m_Vertices) {
    nodeList.add(i->getNodeSetNodeListMeshId(meshID));
  }

  if (order > 1) {
    for (auto i : this->m_Edges) {
      nodeList.add(i->getNodeSetNodeListMeshId(meshID));
    }

    nodeList.add(this->getNodeSetNodeListMeshId(meshID));
  }
  return nodeList;
}

auto LinearQuadrilateralRuntime::getHDivNodes(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  return {};
}

auto LinearQuadrilateralRuntime::getHDivNodesInternal(
    indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  if (order > 0) {
    indexType totnodes = (order - 1) * (order - 1);
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
  return {};
}


auto LinearQuadrilateralRuntime::getH1Shapes(indexType order,
                                             IntegrationPoint &IntegrationPt) -> H1Shapes {

  indexType vertshapes = 4;
  indexType edgeShapes = (order - 1) * 4;
  indexType faceShapes = (order - 1) * (order - 1);

  indexType numshapes = vertshapes + edgeShapes + faceShapes;
  H1Shapes shapes(numshapes, 2);

  prec s1;
  prec ds1;
  prec s2;
  prec ds2;
  prec s3;
  prec ds3;
  prec s4;
  prec ds4;

  //  LobattoShape(s1, ds1, xsi, 0);
  //  LobattoShape(s2, ds2, xsi, 1);
  //  LobattoShape(s3, ds3, eta, 0);
  //  LobattoShape(s4, ds4, eta, 1);
  prec xsi = IntegrationPt.xi;
  prec eta = IntegrationPt.eta;

  LobattoShapes::getShape(s1, ds1, xsi, 0);
  LobattoShapes::getShape(s2, ds2, xsi, 1);
  LobattoShapes::getShape(s3, ds3, eta, 0);
  LobattoShapes::getShape(s4, ds4, eta, 1);

  shapes.shapes(0) = s1 * s3;
  shapes.shapes(1) = s2 * s3;
  shapes.shapes(2) = s2 * s4;
  shapes.shapes(3) = s1 * s4;

  shapes.shapeDeriv(0, 0) = ds1 * s3;
  shapes.shapeDeriv(0, 1) = ds2 * s3;
  shapes.shapeDeriv(0, 2) = ds2 * s4;
  shapes.shapeDeriv(0, 3) = ds1 * s4;

  shapes.shapeDeriv(1, 0) = s1 * ds3;
  shapes.shapeDeriv(1, 1) = s2 * ds3;
  shapes.shapeDeriv(1, 2) = s2 * ds4;
  shapes.shapeDeriv(1, 3) = s1 * ds4;

  if (order > 1) {
    // Edge 1
    indexType counter = 4;
    {
      prec orientation =
          m_Edges[0]->getEdgeOrientation(m_Vertices[0]->getId(), m_Vertices[1]->getId());

      IntegrationPoint eInt;
      eInt.xi = xsi;
      auto eShape = m_Edges[0]->getH1Edge()->getH1ShapesInternal(order, eInt);
      // eShape.shapeDeriv *= orientation;

      for (auto i = 0; i < eShape.shapes.rows(); ++i) {
        if ((i & 1)) {
          eShape.shapeDeriv(0, i) *= orientation;
          eShape.shapes(i) *= orientation;
        }
        shapes.shapes(counter) = eShape.shapes(i) * s3;
        shapes.shapeDeriv(0, counter) = eShape.shapeDeriv(0, i) * s3;
        shapes.shapeDeriv(1, counter) = eShape.shapes(i) * ds3;
        ++counter;
      }
    }
    // Edge 2
    {
      prec orientation = m_Edges[1]->getEdgeOrientation(m_Vertices[1]->getId(),
                                                       m_Vertices[2]->getId());

      IntegrationPoint eInt;
      eInt.xi = eta;
      auto eShape = m_Edges[1]->getH1Edge()->getH1ShapesInternal(order, eInt);
      // eShape.shapeDeriv *= orientation;

      for (auto i = 0; i < eShape.shapes.rows(); ++i) {
        if ((i & 1)) {
          eShape.shapeDeriv(0, i) *= orientation;
          eShape.shapes(i) *= orientation;
        }
        shapes.shapes(counter) = eShape.shapes(i) * s2;
        shapes.shapeDeriv(1, counter) = eShape.shapeDeriv(0, i) * s2;
        shapes.shapeDeriv(0, counter) = eShape.shapes(i) * ds2;
        ++counter;
      }
    }
    // Edge 3
    {
      prec orientation = m_Edges[2]->getEdgeOrientation(m_Vertices[3]->getId(),
                                                       m_Vertices[2]->getId());

      IntegrationPoint eInt;
      eInt.xi = xsi;
      auto eShape = m_Edges[2]->getH1Edge()->getH1ShapesInternal(order, eInt);
      // eShape.shapeDeriv *= orientation;

      for (auto i = 0; i < eShape.shapes.rows(); ++i) {
        if ((i & 1)) {
          eShape.shapeDeriv(0, i) *= orientation;
          eShape.shapes(i) *= orientation;
        }
        shapes.shapes(counter) = eShape.shapes(i) * s4;
        shapes.shapeDeriv(0, counter) = eShape.shapeDeriv(0, i) * s4;
        shapes.shapeDeriv(1, counter) = eShape.shapes(i) * ds4;
        ++counter;
      }
    }
    // Edge 4
    {
      prec orientation = m_Edges[3]->getEdgeOrientation(m_Vertices[0]->getId(),
                                                       m_Vertices[3]->getId());

      IntegrationPoint eInt;
      eInt.xi = eta;
      auto eShape = m_Edges[3]->getH1Edge()->getH1ShapesInternal(order, eInt);
      // eShape.shapeDeriv *= orientation;

      for (auto i = 0; i < eShape.shapes.rows(); ++i) {
        if ((i & 1)) {
          eShape.shapeDeriv(0, i) *= orientation;
          eShape.shapes(i) *= orientation;
        }
        shapes.shapes(counter) = eShape.shapes(i) * s1;
        shapes.shapeDeriv(1, counter) = eShape.shapeDeriv(0, i) * s1;
        shapes.shapeDeriv(0, counter) = eShape.shapes(i) * ds1;
        ++counter;
      }
    }

    auto inShapes = this->getH1ShapesInternal(order, IntegrationPt);
    indexType nn = inShapes.shapes.rows();
    shapes.shapes.block(counter, 0, nn, 1) = inShapes.shapes;
    shapes.shapeDeriv.block(0, counter, 2, nn) = inShapes.shapeDeriv;
    // for (auto i = 0; i < inShapes.shapes.rows(); ++i) {
    //   shapes.shapes(counter) = inShapes.shapes(i);
    //   shapes.shapeDeriv(0, counter) = inShapes.shapeDeriv(0, i);
    //   shapes.shapeDeriv(1, counter) = inShapes.shapeDeriv(1, i);
    //   ++counter;
    // }
  }

  return shapes;
}

auto LinearQuadrilateralRuntime::getH1ShapesInternal(
    indexType order,
    IntegrationPoint &IntegrationPt, faceorientation orientation /*= faceorientation::p_1*/) -> H1Shapes {

  prec xsi = IntegrationPt.xi;
  prec eta = IntegrationPt.eta;

  indexType numShapes = order - 1;
  numShapes *= numShapes;
  H1Shapes shapes(numShapes, 2);
  indexType counter = 0;

  Types::Matrix22<prec> R;

  if (orientation != faceorientation::p_1) {
    R.setZero();
    prec one = prec(1.0);
    switch (orientation) {
    case faceorientation::p_2:
      xsi = -IntegrationPt.eta;
      eta = IntegrationPt.xi;
      R(0, 1) = one;
      R(1, 0) = -one;
      break;
    case faceorientation::p_3:
      xsi = -IntegrationPt.xi;
      eta = -IntegrationPt.eta;
      R(0, 0) = -one;
      R(1, 1) = -one;

      break;
    case faceorientation::p_4:
      xsi = IntegrationPt.eta;
      eta = -IntegrationPt.xi;
      R(0, 1) = -one;
      R(1, 0) = one;
      break;
    case faceorientation::n_1:
      xsi = -IntegrationPt.xi;
      eta = IntegrationPt.eta;
      R(0, 0) = -one;
      R(1, 1) = one;
      break;
    case faceorientation::n_2:
      xsi = -IntegrationPt.eta;
      eta = -IntegrationPt.xi;
      R(0, 1) = -one;
      R(1, 0) = -one;
      break;
    case faceorientation::n_3:
      xsi = IntegrationPt.xi;
      eta = -IntegrationPt.eta;
      R(0, 0) = one;
      R(1, 1) = -one;
      break;
    case faceorientation::n_4:
      xsi = IntegrationPt.eta;
      eta = IntegrationPt.xi;
      R(0, 1) = one;
      R(1, 0) = one;
      break;
    default: {
    }
    }
  }

  for (auto i = 2; i <= order; i++) {
    auto xiShape = LobattoShapes::getShape(xsi, i);
    // xiShape.shapeDerivative *= factXi;
    for (auto j = 2; j <= order; j++) {
      auto etaShape = LobattoShapes::getShape(eta, j);
      // etaShape.shapeDerivative *= factEta;
      shapes.shapes(counter) = xiShape.shapeValue * etaShape.shapeValue;
      shapes.shapeDeriv(0, counter) =
          xiShape.shapeDerivative * etaShape.shapeValue;
      shapes.shapeDeriv(1, counter) =
          xiShape.shapeValue * etaShape.shapeDerivative;
      ++counter;
    }
  }
  if (orientation != faceorientation::p_1) {
    shapes.shapeDeriv = R * shapes.shapeDeriv;
  }

  return shapes;
}

// HDivShapes
void LinearQuadrilateralRuntime::setHDivShapes(indexType meshId,
                                               indexType order,
                                               NodeTypes type) {
  indexType numNodes = order + 1;
  for (auto i = 0; i < 4; ++i) {
    m_Edges[i]->setNodeSet(meshId, numNodes, type);
  }

  if (order >= 1) {
    numNodes = 2 * order * order + 2 * order;
    this->setNodeSet(meshId, numNodes, type);
  }
}

void LinearQuadrilateralRuntime::getHDivDofs(
    std::vector<DegreeOfFreedom *> &Dofs,
    indexType meshID, indexType order, NodeTypes type) {
  for (auto i = 0; i < 4; ++i) {
    auto nodeList = m_Edges[i]->getNodeSetNodeListMeshId(meshID);
    auto tdofs(nodeList.getDegreesOfFreedom());
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }

  if (order >= 1) {
    auto nodeList = this->getNodeSetNodeListMeshId(meshID);
    auto tdofs(nodeList.getDegreesOfFreedom());
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
}

void LinearQuadrilateralRuntime::getHDivShapes(indexType order,
                                               Types::Matrix2X<prec> &shape,
                                               Types::VectorX<prec> &dshape,
                                               prec xi, prec eta) {
  indexType edgeShapes = (order + 1) * 4;
  indexType faceShapes = 2 * order + 2 * order * order;

  indexType numshapes = edgeShapes + faceShapes;
  shape.resize(2, numshapes);
  dshape.resize(numshapes);

  prec l0x;
  prec dl0x;
  prec l1x;
  prec dl1x;
  prec l0e;
  prec dl0e;
  prec l1e;
  prec dl1e;
  prec Lkx;
  prec dLkx;
  prec Lke;
  prec dLke;
  LobattoShapes::getShape(l0x, dl0x, xi, 0);
  LobattoShapes::getShape(l1x, dl1x, xi, 1);
  LobattoShapes::getShape(l0e, dl0e, eta, 0);
  LobattoShapes::getShape(l1e, dl1e, eta, 1);

  int counter = 0;

  // edge1
  {
    prec orientation = m_Edges[0]->getEdgeOrientation(m_Vertices[0]->getId(),
                                                     m_Vertices[1]->getId());
    prec orientationx = orientation;
    for (auto i = 0; i <= order; i++) {
      LegendreShapes::getShape(Lkx, dLkx, orientationx * xi, i);

      shape(0, counter) = 0;
      shape(1, counter) = -1 * l0e * Lkx * orientation;
      dshape(counter) = -1 * dl0e * Lkx * orientation; // * orientationx;
      ++counter;
    }
  }

  // edge2
  {
    prec orientation = m_Edges[1]->getEdgeOrientation(m_Vertices[1]->getId(),
                                                     m_Vertices[2]->getId());
    prec orientationx = orientation;
    for (auto i = 0; i <= order; i++) {
      LegendreShapes::getShape(Lke, dLke, orientationx * eta, i);

      shape(0, counter) = l1x * Lke * orientation;
      shape(1, counter) = 0;
      dshape(counter) = dl1x * Lke * orientation; // * orientationx;
      ++counter;
    }
  }

  // edge3
  {
    prec orientation = m_Edges[2]->getEdgeOrientation(m_Vertices[2]->getId(),
                                                     m_Vertices[3]->getId());
    prec orientationx = orientation;
    for (auto i = 0; i <= order; i++) {
      LegendreShapes::getShape(Lkx, dLkx, orientationx * xi, i);

      shape(0, counter) = 0;
      shape(1, counter) = l1e * Lkx * orientation;
      dshape(counter) = dl1e * Lkx * orientation; // * orientationx;
      ++counter;
    }
  }

  // edge4
  {
    prec orientation = m_Edges[3]->getEdgeOrientation(m_Vertices[3]->getId(),
                                                     m_Vertices[0]->getId());
    prec orientationx = orientation;
    for (auto i = 0; i <= order; i++) {
      LegendreShapes::getShape(Lke, dLke, orientationx * eta, i);

      shape(0, counter) = -1 * l0x * Lke * orientation;
      shape(1, counter) = 0;
      dshape(counter) = -1 * dl0x * Lke * orientation; // *orientationx;
      ++counter;
    }
  }

  // Face
  prec lkx;
  prec dlkx;
  prec lke;
  prec dlke;

  if (order >= 1) {
    for (auto a = 2; a <= order + 1; a++) {
      for (auto b = 0; b <= order; b++) {
        LobattoShapes::getShape(lkx, dlkx, xi, a);
        LegendreShapes::getShape(Lke, dLke, eta, b);
        LobattoShapes::getShape(lke, dlke, eta, a);
        LegendreShapes::getShape(Lkx, dLkx, xi, b);

        shape(0, counter) = lkx * Lke;
        shape(1, counter) = 0;
        dshape(counter) = dlkx * Lke;
        ++counter;
        shape(0, counter) = 0;
        shape(1, counter) = Lkx * lke;
        dshape(counter) = Lkx * dlke;
        ++counter;
      }
    }
  }
}

auto LinearQuadrilateralRuntime::getHDivShapes(indexType order,
                                               IntegrationPoint &IntegrationPt)
    -> HDivShapes {
  HDivShapes shapes;

  prec xi = IntegrationPt.xi;
  prec eta = IntegrationPt.eta;

  indexType edgeShapes = (order + 1) * 4;
  indexType faceShapes = 2 * order + 2 * order * order;

  indexType numshapes = edgeShapes + faceShapes;
  shapes.shapes.resize(2, numshapes);
  shapes.shapeDeriv.resize(1, numshapes);

  prec l0x;
  prec dl0x;
  prec l1x;
  prec dl1x;
  prec l0e;
  prec dl0e;
  prec l1e;
  prec dl1e;
  prec Lkx;
  prec dLkx;
  prec Lke;
  prec dLke;
  LobattoShapes::getShape(l0x, dl0x, xi, 0);
  LobattoShapes::getShape(l1x, dl1x, xi, 1);
  LobattoShapes::getShape(l0e, dl0e, eta, 0);
  LobattoShapes::getShape(l1e, dl1e, eta, 1);

  int counter = 0;

  // edge1
  {
    prec orientation = m_Edges[0]->getEdgeOrientation(m_Vertices[0]->getId(),
                                                     m_Vertices[1]->getId());
    prec orientationx = orientation;
    for (auto i = 0; i <= order; i++) {
      LegendreShapes::getShape(Lkx, dLkx, orientationx * xi, i);
      // if (i % 2 == 0) {
      shapes.shapes(0, counter) = prec(0);
      shapes.shapes(1, counter) = -prec(1) * l0e * Lkx * orientation;
      shapes.shapeDeriv(counter) =
          -1 * dl0e * Lkx *
          orientation; // * orientationx;         //unter Vorbehalt!!
      //}
      // else {
      //    shapes.shapes(0, counter) = 0;
      //    shapes.shapes(1, counter) = -1 * l0e * Lkx ;
      //    shapes.shapeDeriv(counter) = -1 * dl0e * Lkx;         //unter
      //    Vorbehalt!!
      //}
      ++counter;
    }
  }

  // edge2
  {
    prec orientation = m_Edges[1]->getEdgeOrientation(m_Vertices[1]->getId(),
                                                     m_Vertices[2]->getId());
    prec orientationx = orientation;
    for (auto i = 0; i <= order; i++) {
      LegendreShapes::getShape(Lke, dLke, orientationx * eta, i);
      // if (i % 2 == 0) {
      shapes.shapes(0, counter) = l1x * Lke * orientation;
      shapes.shapes(1, counter) = prec(0);
      shapes.shapeDeriv(counter) =
          dl1x * Lke *
          orientation; // * orientationx;              //unter Vorbehalt!!
      //}
      // else {
      //    shapes.shapes(0, counter) = l1x * Lke;
      //    shapes.shapes(1, counter) = 0;
      //    shapes.shapeDeriv(counter) = dl1x * Lke * orientation; //unter
      //    Vorbehalt!!
      //}
      ++counter;
    }
  }

  // edge3
  {
    prec orientation = m_Edges[2]->getEdgeOrientation(m_Vertices[2]->getId(),
                                                     m_Vertices[3]->getId());
    prec orientationx = -orientation;
    for (auto i = 0; i <= order; i++) {
      LegendreShapes::getShape(Lkx, dLkx, orientationx * xi, i);
      // if (i % 2 == 0) {
      shapes.shapes(0, counter) = prec(0);
      shapes.shapes(1, counter) = l1e * Lkx * orientation;
      shapes.shapeDeriv(counter) =
          dl1e * Lkx *
          orientation; // * orientationx;              //unter Vorbehalt!!
      //}
      // else {
      //    shapes.shapes(0, counter) = 0;
      //    shapes.shapes(1, counter) = -1 * l1e * Lkx;
      //    shapes.shapeDeriv(counter) = dl1e * Lkx * orientation; //unter
      //    Vorbehalt!!
      //}
      ++counter;
    }
  }

  // edge4
  {
    prec orientation = m_Edges[3]->getEdgeOrientation(m_Vertices[3]->getId(),
                                                     m_Vertices[0]->getId());
    prec orientationx = -orientation;
    for (auto i = 0; i <= order; i++) {
      LegendreShapes::getShape(Lke, dLke, orientationx * eta, i);
      // if (i % 2 == 0) {
      shapes.shapes(0, counter) = -prec(1) * l0x * Lke * orientation;
      shapes.shapes(1, counter) = prec(0);
      shapes.shapeDeriv(counter) =
          -prec(1) * dl0x * Lke *
          orientation; // *orientationx;                  //unter Vorbehalt!!
      //}
      // else {
      //    shapes.shapes(0, counter) = l0x * Lke;
      //    shapes.shapes(1, counter) = 0;
      //    shapes.shapeDeriv(counter) = -1 * dl0x * Lke * orientation; //unter
      //    Vorbehalt!!
      //}
      ++counter;
    }
  }

  // Face
  prec lkx;
  prec dlkx;
  prec lke;
  prec dlke;

  if (order >= 1) {
    for (auto a = 2; a <= order + 1; a++) {
      for (auto b = 0; b <= order; b++) {
        LobattoShapes::getShape(lkx, dlkx, xi, a);
        LegendreShapes::getShape(Lke, dLke, eta, b);
        LobattoShapes::getShape(lke, dlke, eta, a);
        LegendreShapes::getShape(Lkx, dLkx, xi, b);

        shapes.shapes(0, counter) = lkx * Lke;
        shapes.shapes(1, counter) = 0;
        shapes.shapeDeriv(counter) = dlkx * Lke; // unter Vorbehalt!!
        ++counter;
        shapes.shapes(0, counter) = 0;
        shapes.shapes(1, counter) = Lkx * lke;
        shapes.shapeDeriv(counter) = Lkx * dlke; // unter Vorbehalt!!
        ++counter;
      }
    }
  }

  return shapes;
}

void LinearQuadrilateralRuntime::setAllNodeBoundaryConditionMeshId(
    indexType meshId, indexType dof) {
  GeometryBaseRuntime::setAllNodeBoundaryConditionMeshId(meshId, dof);

  for (auto &i : m_Edges) {
    i->setAllNodeBoundaryConditionMeshId(meshId, dof);
  }
}

void LinearQuadrilateralRuntime::geometryToParaview(
    vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh) {
  indexType numPoints = 4;
  std::vector<indexType> points(numPoints);
  points.clear();
  for (auto &i : m_Vertices) {
    i->geometryToParaview(paraviewAdapter, mainMesh, subMesh);
    points.push_back(i->getId());
  }

  paraviewAdapter.addCell(mainMesh, subMesh, this->m_Face_Data_Element.getId(), 1, points, 4, VTK_QUAD);
}

void LinearQuadrilateralRuntime::computeWeightsParaview(
    vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh) {

  auto GP = this->getIntegrationPoints(-1);
  GP.setOrder(2);

  for (auto i : GP) {
    auto jaco = this->getJacobian(i);
    auto shapes = this->getH1Shapes(1, i);
    prec dA = jaco.determinant() * i.weight;

    for (auto j = 0; j < 4; ++j) {
      std::vector<prec> val;
      val.push_back(shapes.shapes(j) * dA);
      paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, val,
                                           m_Vertices[j]->getId(), 1,
                                           paraviewNames::weightName());
    }
  }
}

void LinearQuadrilateralRuntime::H1SolutionToParaview(
    vtkPlotInterface &paraviewAdapter, indexType mainMesh, indexType subMesh,
    indexType order, Types::VectorX<prec> &solution, std::string &name) {

  for (auto i = 0; i < 4; ++i) {
    std::vector<prec> sol(3);
    for (auto j = 0; j < 3; ++j) {
      sol[j] = solution(3 * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, m_Vertices[i]->getId(),
                                 sol, 3, name);
  }
}

void LinearQuadrilateralRuntime::H1DataToParaview(
    vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, Types::VectorX<prec> &Data,
    indexType numberComponents, indexType order, std::string &name) {
  for (auto i = 0; i < 4; ++i) {
    std::vector<prec> sol(numberComponents);
    for (auto j = 0; j < numberComponents; ++j) {
      sol[j] = Data(numberComponents * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, m_Vertices[i]->getId(),
                                 sol, numberComponents, name);
  }
}

void LinearQuadrilateralRuntime::projectDataToParaviewVertices(
    vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {

  auto shapes = this->getH1Shapes(1, IntegrationPt);
  std::vector<prec> vals(numberComponents);

  auto jaco = this->getJacobian(IntegrationPt);
  auto dA = jaco.determinant() * IntegrationPt.weight;
  for (auto i = 0; i < 4; ++i) {
    for (auto j = 0; j < numberComponents; ++j) {
      vals[j] = data(j) * shapes.shapes(i) * dA;
    }
    paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, vals,
                                         m_Vertices[i]->getId(),
                                         numberComponents, name);
  }
}

void LinearQuadrilateralRuntime::setLoad(
    LoadList &loadlist, indexType meshid, ShapeFunctionTypes shapeType,
    indexType shapeOrder, Types::VectorX<prec> &Loads, indexType propNumber,
    Types::VectorX<prec> &direction, bool local, bool add) {

  if (shapeType == ShapeFunctionTypes::H1) {
    auto GP = this->getIntegrationPoints(-1);
    GP.setOrder(shapeOrder * 2);
    auto nodes = this->getH1Nodes(meshid, shapeOrder);
    Types::VectorX<prec> tempLoads(nodes.size() * 3);
    tempLoads.setZero();
    for (auto i : GP) {
      auto shape = this->getH1Shapes(shapeOrder, i);
      Types::Vector3<prec> dx;
      Types::Vector3<prec> dy;
      dx.setZero();
      dy.setZero();
      for (auto j = 0; j < 4; ++j) {
        dx += shape.shapeDeriv(0, j) * m_Vertices[j]->getCoordinates();
        dy += shape.shapeDeriv(1, j) * m_Vertices[j]->getCoordinates();
      }
      auto dA = dx.cross(dy).norm() * i.weight;
      indexType counter = 0;
      for (auto j = 0; j < nodes.size(); ++j) {
        auto tempDofs = nodes[j]->getDegreesOfFreedom();
        auto daq = shape.shapes(j) * dA;
        for (auto k = 0; k < tempDofs.size(); ++k) {
          auto ll = daq * Loads(k);
          tempLoads(counter) += ll;
          counter++;
        }
      }
    }
    indexType counter = 0;
    for (auto i : nodes) {
      auto tempDofs = i->getDegreesOfFreedom();
      for (auto j : tempDofs) {
        loadlist.setLoad(propNumber, j->getId(),
                                        tempLoads(counter), add);
        counter++;
      }
    }
  }
}

void LinearQuadrilateralRuntime::setPrescribedSolution(
    LoadList &loadlist, indexType meshid, ShapeFunctionTypes shapeType,
    indexType shapeOrder, Types::VectorX<prec> &Solution, indexType propNumber,
    Types::VectorX<prec> &direction, bool local, bool add) {
  if (shapeType == ShapeFunctionTypes::H1) {
    for (auto &i : m_Vertices) {
      i->setPrescribedSolution(loadlist, meshid, shapeType, shapeOrder,
                               Solution, propNumber, direction, local, add);
    }
  }
}

void LinearQuadrilateralRuntime::setBoundaryCondition(
    indexType meshId, indexType order,
    ShapeFunctionTypes shapeType, Types::Vector3<indexType> &dofs, bool set) {
  if (shapeType == ShapeFunctionTypes::H1) {
    auto nodes = this->getH1Nodes(meshId, order);
    for (auto &i : nodes) {
      if (set) {
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            i->setBoundaryCondition(j);
          } else {
            i->unsetBoundaryCondition(j);
          }
        }
      } else {
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            i->setBoundaryCondition(j);
          }
        }
      }
    }
  }
}

void LinearQuadrilateralRuntime::setSpecialPlateShapes(
    indexType meshid, indexType order,
    NodeTypes type) {}
auto LinearQuadrilateralRuntime::getSpecialPlateDofs(
    indexType meshID, indexType order,
    NodeTypes type) -> std::vector<DegreeOfFreedom *> {
  std::vector<DegreeOfFreedom *> Dofs;
  // auto &edge = pointers.getGeometryData()->getEdge(this->edges[0]);
  return Dofs;
}
auto LinearQuadrilateralRuntime::getSpecialPlateShapes(
    IntegrationPoint &intPoint, indexType order)
    -> SpecialPlateShapes {
  SpecialPlateShapes shapes;
  prec fValue;
  prec fDeriv;
  LobattoShapes::getShape(fValue, fDeriv, intPoint.xi, 0);
  return shapes;
}

void LinearQuadrilateralRuntime::checkUpdateElement(EquationHandler &eqHandler,
                                                    GeometryData &geoData) {
  m_Face_Data_Element.checkUpdateElement(eqHandler,geoData);
}


void LinearQuadrilateralRuntime::set_geometry_pointers(GeometryData &geoData) {
  m_Face_Data_Element.set_geometry_pointers(geoData);
}



const GeometryTypes LinearQuadrilateralRuntime::type =
    GeometryTypes::LinearQuadrilateral;
} // namespace HierAMuS::Geometry
