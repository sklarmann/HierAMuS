// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include "geometry/GeometryData.h"
#include "MatrixTypes.h"
#include "geometry/Faces/FacesData.h"
#include "geometry/Faces/LinearQuadrilateralRuntime.h"
#include "geometry/GeometryBaseData.h"
#include "geometry/GeometryTypes.h"
#include "plot/vtkplotClass.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "shapefunctions/LegendreShapes.h"
#include "shapefunctions/LobattoShapes.h"
#include <exception>


#include "GenericNodes.h"

#include <geometry/Faces/LinearQuadrilateralData.h>
#include <sstream>

#include <geometry/Edges/EdgesData.h>
#include <geometry/VertexData.h>

#include "LoadList.h"

#include <stdexcept>
#include <vector>
#include <vtkCellType.h>

#include "HelperFunctions.h"

#include "shapefunctions/IntegrationsPoints/helperClasses/IntegrationPointsManagement.h"

namespace HierAMuS::Geometry {

LinearQuadrilateralData::LinearQuadrilateralData() : FacesDataInterface(){};

LinearQuadrilateralData::~LinearQuadrilateralData() = default;

auto LinearQuadrilateralData::getRuntimeObject(GeometryData &geoData)
    -> std::shared_ptr<FacesRuntime> {
  return std::make_shared<LinearQuadrilateralRuntime>(geoData, *this);
}

auto LinearQuadrilateralData::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::LinearQuadrilateralData::type;
}

auto LinearQuadrilateralData::getCoordinates(prec xi, prec eta)
    -> Types::Vector3<prec> {
  Types::Vector3<prec> ret;
  Types::Vector3<prec> t1;
  ret.setZero();
  prec s1;
  prec ds1;
  prec s2;
  prec ds2;
  prec s3;
  prec ds3;
  prec s4;
  prec ds4;
  std::array<prec, 4> vals;
  //  LobattoShape(s1, ds1, xi, 0);
  //  LobattoShape(s2, ds2, xi, 1);
  //  LobattoShape(s3, ds3, eta, 0);
  //  LobattoShape(s4, ds4, eta, 1);

  LobattoShapes::getShape(s1, ds1, xi, 0);
  LobattoShapes::getShape(s2, ds2, xi, 1);
  LobattoShapes::getShape(s3, ds3, eta, 0);
  LobattoShapes::getShape(s4, ds4, eta, 1);

  vals[0] = s1 * s3;
  vals[1] = s2 * s3;
  vals[2] = s2 * s4;
  vals[3] = s1 * s4;

  for (auto i = 0; i < 4; ++i) {
    t1 = m_verts_pointers[i]->getCoordinates();
    ret += t1 * vals[i];
  }
  return ret;
}

auto LinearQuadrilateralData::getCoordinates(IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  Types::Vector3<prec> ret;
  Types::Vector3<prec> t1;
  ret.setZero();
  prec s1;
  prec ds1;
  prec s2;
  prec ds2;
  prec s3;
  prec ds3;
  prec s4;
  prec ds4;
  std::array<prec, 4> vals;
  //  LobattoShape(s1, ds1, xi, 0);
  //  LobattoShape(s2, ds2, xi, 1);
  //  LobattoShape(s3, ds3, eta, 0);
  //  LobattoShape(s4, ds4, eta, 1);

  LobattoShapes::getShape(s1, ds1, integrationPoint.xi, 0);
  LobattoShapes::getShape(s2, ds2, integrationPoint.xi, 1);
  LobattoShapes::getShape(s3, ds3, integrationPoint.eta, 0);
  LobattoShapes::getShape(s4, ds4, integrationPoint.eta, 1);

  vals[0] = s1 * s3;
  vals[1] = s2 * s3;
  vals[2] = s2 * s4;
  vals[3] = s1 * s4;

  for (auto i = 0; i < 4; ++i) {
    t1 = m_verts_pointers[i]->getCoordinates();
    ret += t1 * vals[i];
  }
  return ret;
}

auto LinearQuadrilateralData::getTangent_G1(IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  // auto shapes = this->getH1Shapes(pointers, 1, integrationPoint);
  Types::Vector3<prec> ret;

  auto [x1, dx1] = LobattoShapes::getShape(integrationPoint.xi, 0);
  auto [y1, dy1] = LobattoShapes::getShape(integrationPoint.eta, 0);
  auto [x2, dx2] = LobattoShapes::getShape(integrationPoint.xi, 1);
  auto [y2, dy2] = LobattoShapes::getShape(integrationPoint.eta, 1);

  ret = m_verts_pointers[0]->getCoordinates() * (dx1 * y1);
  ret += m_verts_pointers[1]->getCoordinates() * (dx2 * y1);
  ret += m_verts_pointers[2]->getCoordinates() * (dx2 * y2);
  ret += m_verts_pointers[3]->getCoordinates() * (dx1 * y2);

  return ret;
}

auto LinearQuadrilateralData::getTangent_G2(IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  // auto shapes = this->getH1Shapes(pointers, 1, integrationPoint);
  Types::Vector3<prec> ret;

  auto [x1, dx1] = LobattoShapes::getShape(integrationPoint.xi, 0);
  auto [y1, dy1] = LobattoShapes::getShape(integrationPoint.eta, 0);
  auto [x2, dx2] = LobattoShapes::getShape(integrationPoint.xi, 1);
  auto [y2, dy2] = LobattoShapes::getShape(integrationPoint.eta, 1);

  ret = m_verts_pointers[0]->getCoordinates() * (x1 * dy1);
  ret += m_verts_pointers[1]->getCoordinates() * (x2 * dy1);
  ret += m_verts_pointers[2]->getCoordinates() * (x2 * dy2);
  ret += m_verts_pointers[3]->getCoordinates() * (x1 * dy2);
  return ret;
}

auto LinearQuadrilateralData::getFaceNormal()
    -> Types::Vector3<prec> {
  auto &V1 = *m_verts_pointers[0];
  auto &V2 = *m_verts_pointers[1];
  auto &V3 = *m_verts_pointers[3];

  Types::Vector3<prec> dx = V2.getCoordinates() - V1.getCoordinates();
  Types::Vector3<prec> dy = V3.getCoordinates() - V1.getCoordinates();

  Types::Vector3<prec> n = dx.cross(dy).normalized();
  return n;
}

auto LinearQuadrilateralData::getOrientation(indexType vertex1,
                                             indexType vertex2)
    -> faceorientation {
  faceorientation result = faceorientation::p_1;
  // look for first vertex position
  indexType i = 0;
  if (this->m_verts[0] == vertex1) {
    i = 0;
  } else if (this->m_verts[1] == vertex1) {
    i = 1;
  } else if (this->m_verts[2] == vertex1) {
    i = 2;
  } else if (this->m_verts[3] == vertex1) {
    i = 3;
  } else {
    throw std::runtime_error(
        "Vertex1 not found in faceorientation of linearquadrilateral");
  }

  // look for second vertex position
  indexType nextVertex = i + 1;
  indexType prevVertex = i - 1;
  if (nextVertex > 3)
    nextVertex = 0;
  if (prevVertex < 0)
    prevVertex = 3;
  if (this->m_verts[nextVertex] == vertex2) {
    if (i == 0) {
      result = faceorientation::p_1;
    } else if (i == 1) {
      result = faceorientation::p_2;
    } else if (i == 2) {
      result = faceorientation::p_3;
    } else if (i == 3) {
      result = faceorientation::p_4;
    }
  } else if (this->m_verts[prevVertex] == vertex2) {
    if (i == 1) {
      result = faceorientation::n_1;
    } else if (i == 2) {
      result = faceorientation::n_2;
    } else if (i == 3) {
      result = faceorientation::n_3;
    } else if (i == 0) {
      result = faceorientation::n_4;
    }
  } else {
    throw std::runtime_error(
        "Something wrong in faceorientation for linearquadrilateral");
  }
  // result = faceorientation::p_1;
  return result;
}

void LinearQuadrilateralData::modifyIntegrationpoint(
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

auto LinearQuadrilateralData::getIntegrationPoints(indexType elementId)
    -> IntegrationPoints {
  IntegrationPoints temp = IntegrationPointsManagement::getIntegrationsPoints(elementId);
  temp.setType(IntegrationType::Gauss2D);
  return temp;
}

auto LinearQuadrilateralData::getJacobian(IntegrationPoint &IntegrationPt)
    -> Types::Matrix22<prec> {
  Types::Matrix22<prec> jacobi;

  prec N1x = IntegrationPt.eta * prec(0.25) - prec(0.25);
  prec N2x = prec(0.25) - IntegrationPt.eta * prec(0.25);
  prec N3x = IntegrationPt.eta * prec(0.25) + prec(0.25);
  prec N4x = -IntegrationPt.eta * prec(0.25) - prec(0.25);

  prec N1e = IntegrationPt.xi * prec(0.25) - prec(0.25);
  prec N2e = -IntegrationPt.xi * prec(0.25) - prec(0.25);
  prec N3e = IntegrationPt.xi * prec(0.25) + prec(0.25);
  prec N4e = prec(0.25) - IntegrationPt.xi * prec(0.25);

  // Vertex 1
  auto coord = m_verts_pointers[0]->getCoordinates();
  jacobi(0, 0) = N1x * coord(0);
  jacobi(1, 0) = N1x * coord(1);
  jacobi(0, 1) = N1e * coord(0);
  jacobi(1, 1) = N1e * coord(1);

  // Vertex 2
  coord = m_verts_pointers[1]->getCoordinates();
  jacobi(0, 0) += N2x * coord(0);
  jacobi(1, 0) += N2x * coord(1);
  jacobi(0, 1) += N2e * coord(0);
  jacobi(1, 1) += N2e * coord(1);

  // Vertex 3
  coord = m_verts_pointers[2]->getCoordinates();
  jacobi(0, 0) += N3x * coord(0);
  jacobi(1, 0) += N3x * coord(1);
  jacobi(0, 1) += N3e * coord(0);
  jacobi(1, 1) += N3e * coord(1);

  // Vertex 4
  coord = m_verts_pointers[3]->getCoordinates();
  jacobi(0, 0) += N4x * coord(0);
  jacobi(1, 0) += N4x * coord(1);
  jacobi(0, 1) += N4e * coord(0);
  jacobi(1, 1) += N4e * coord(1);

  return jacobi;
}

// H1Shapes

void LinearQuadrilateralData::setH1Shapes(indexType meshId, indexType order,
                                          NodeTypes type) {
  for (auto i = 0; i < 4; ++i) {
    m_edges_pointers[i]->setH1Shapes(meshId, order, type);
  }
  this->setH1ShapesInternal(meshId, order, type);
}

void LinearQuadrilateralData::setH1ShapesInternal(indexType meshId,
                                                  indexType order,
                                                  NodeTypes type) {
  if (order > 1) {
    indexType numNodes = order - 1;
    numNodes *= numNodes;

    this->setNodeSet(meshId, numNodes, type);
  }
}

auto LinearQuadrilateralData::getH1Dofs(indexType meshID, indexType order)
    -> std::vector<DegreeOfFreedom *> {
  // std::vector<DegreeOfFreedom *> Dofs;
  // this->getH1Dofs(pointers, Dofs, meshID, order);
  // return Dofs;
  auto nodeList = this->getH1NodesList(meshID, order);
  return nodeList.getDegreesOfFreedom();
}

void LinearQuadrilateralData::getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs,
                                        indexType meshID, indexType order) {
  for (auto i = 0; i < 4; ++i) {
    auto nodeList = m_verts_pointers[i]->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }
  if (order > 1) {
    for (auto i = 0; i < 4; ++i) {
      m_edges_pointers[i]->getH1DofsInternal(Dofs, meshID, order);
    }
    this->getH1DofsInternal(Dofs, meshID, order);
  }
}

void LinearQuadrilateralData::getH1DofsInternal(
    std::vector<DegreeOfFreedom *> &Dofs,
    indexType meshID, indexType order) {
  if (order > 1) {
    auto nodeList = this->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }
}

auto LinearQuadrilateralData::getH1Nodes(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  for (auto i = 0; i < 4; ++i) {
    auto tempNodes = m_verts_pointers[i]->getNodesOfSet(meshID);
    nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
  }
  if (order > 1) {
    for (auto i = 0; i < 4; ++i) {
      auto tempNodes =
          m_edges_pointers[i]->getH1NodesInternal(meshID, order);
      nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
    }
    auto tempNodes = this->getH1NodesInternal(meshID, order);
    nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
  }
  return nodes;
}

auto LinearQuadrilateralData::getH1NodesInternal(indexType meshID,
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

auto LinearQuadrilateralData::getH1NodesList(indexType meshID, indexType order)
    -> MeshIdNodeList {

  MeshIdNodeList nodeList(meshID);
  nodeList.reserve(9);
  for (auto i : this->m_verts_pointers) {
    auto templist = i->getNodeSetNodeListMeshId(meshID);
    nodeList.add(templist);
  }
  if (order > 1) {
    for (auto i : this->m_edges_pointers) {
      auto templist = i->getNodeSetNodeListMeshId(meshID);
      nodeList.add(templist);
    }

    auto tempList = this->getNodeSetNodeListMeshId(meshID);
    nodeList.add(tempList);
  }
  return nodeList;
}

auto LinearQuadrilateralData::getHDivNodes(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  return {};
}

auto LinearQuadrilateralData::getHDivNodesInternal(indexType meshID,
                                                   indexType order)
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

void LinearQuadrilateralData::getH1Shapes(
    indexType order, Types::VectorX<prec> &shape,
    Types::Matrix2X<prec> &shapeDerivative, prec xsi, prec eta) {

  indexType vertshapes = 4;
  indexType edgeShapes = (order - 1) * 4;
  indexType faceShapes = (order - 1) * (order - 1);

  indexType numshapes = vertshapes + edgeShapes + faceShapes;
  shape.resize(numshapes);
  shapeDerivative.resize(2, numshapes);

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

  LobattoShapes::getShape(s1, ds1, xsi, 0);
  LobattoShapes::getShape(s2, ds2, xsi, 1);
  LobattoShapes::getShape(s3, ds3, eta, 0);
  LobattoShapes::getShape(s4, ds4, eta, 1);

  shape(0) = s1 * s3;
  shape(1) = s2 * s3;
  shape(2) = s2 * s4;
  shape(3) = s1 * s4;

  shapeDerivative(0, 0) = ds1 * s3;
  shapeDerivative(0, 1) = ds2 * s3;
  shapeDerivative(0, 2) = ds2 * s4;
  shapeDerivative(0, 3) = ds1 * s4;

  shapeDerivative(1, 0) = s1 * ds3;
  shapeDerivative(1, 1) = s2 * ds3;
  shapeDerivative(1, 2) = s2 * ds4;
  shapeDerivative(1, 3) = s1 * ds4;

  if (order > 1) {
    Types::VectorX<prec> tempshape;
    Types::VectorX<prec> tempEdgeDerivative;
    Types::Matrix2X<prec> tempshapeDerivative;
    // Edge 1
    indexType counter = 4;
    {
      auto &tempEdge = *m_edges_pointers[0];
      prec orientation =
          tempEdge.getEdgeOrientation(this->m_verts[0], this->m_verts[1]);
      // prec txsi = xsi * orientation;
      tempEdge.getH1ShapesInternal(order, tempshape,
                                   tempEdgeDerivative, xsi);
      // tempEdgeDerivative *= orientation;
      for (auto i = 0; i < tempshape.rows(); ++i) {
        if ((i & 1)) {
          tempEdgeDerivative(i) *= orientation;
          tempshape(i) *= orientation;
        }
        shape(counter) = tempshape(i) * s3;
        shapeDerivative(0, counter) = tempEdgeDerivative(i) * s3;
        shapeDerivative(1, counter) = tempshape(i) * ds3;
        ++counter;
      }
    }
    // Edge 2
    {
      auto &tempEdge = *m_edges_pointers[1];
      prec orientation =
          tempEdge.getEdgeOrientation(this->m_verts[1], this->m_verts[2]);
      // prec txsi = eta * orientation;
      tempEdge.getH1ShapesInternal(order, tempshape,
                                   tempEdgeDerivative, eta);
      // tempEdgeDerivative *= orientation;
      for (auto i = 0; i < tempshape.rows(); ++i) {
        if ((i & 1)) {
          tempEdgeDerivative(i) *= orientation;
          tempshape(i) *= orientation;
        }
        shape(counter) = tempshape(i) * s2;
        shapeDerivative(1, counter) = tempEdgeDerivative(i) * s2;
        shapeDerivative(0, counter) = tempshape(i) * ds2;
        ++counter;
      }
    }
    // Edge 3
    {
      auto &tempEdge = *m_edges_pointers[2];
      prec orientation =
          tempEdge.getEdgeOrientation(this->m_verts[3], this->m_verts[2]);
      // prec txsi = xsi * orientation;
      tempEdge.getH1ShapesInternal(order, tempshape,
                                   tempEdgeDerivative, xsi);
      // tempEdgeDerivative *= orientation;
      for (auto i = 0; i < tempshape.rows(); ++i) {
        if ((i & 1)) {
          tempEdgeDerivative(i) *= orientation;
          tempshape(i) *= orientation;
        }
        shape(counter) = tempshape(i) * s4;
        shapeDerivative(0, counter) = tempEdgeDerivative(i) * s4;
        shapeDerivative(1, counter) = tempshape(i) * ds4;
        ++counter;
      }
    }
    // Edge 4
    {
      auto &tempEdge = *m_edges_pointers[3];
      prec orientation =
          tempEdge.getEdgeOrientation(this->m_verts[0], this->m_verts[3]);
      // prec txsi = eta * orientation;
      tempEdge.getH1ShapesInternal(order, tempshape,
                                   tempEdgeDerivative, eta);
      // tempEdgeDerivative *= orientation;
      for (auto i = 0; i < tempshape.rows(); ++i) {
        if ((i & 1)) {
          tempEdgeDerivative(i) *= orientation;
          tempshape(i) *= orientation;
        }
        shape(counter) = tempshape(i) * s1;
        shapeDerivative(1, counter) = tempEdgeDerivative(i) * s1;
        shapeDerivative(0, counter) = tempshape(i) * ds1;
        ++counter;
      }
    }
    this->getH1ShapesInternal(order, tempshape, tempshapeDerivative,
                              xsi, eta);
    indexType nn = tempshape.cols();
    std::cout << nn << std::endl;
    for (auto i = 0; i < tempshape.cols(); ++i) {
      shape(counter) = tempshape(i);
      shapeDerivative(0, counter) = tempshapeDerivative(0, i);
      shapeDerivative(1, counter) = tempshapeDerivative(1, i);
      ++counter;
    }
  }
}

void LinearQuadrilateralData::getH1ShapesInternal(
    indexType order, Types::VectorX<prec> &shape,
    Types::Matrix2X<prec> &shapeDerivative, prec xsi, prec eta) {

  if (order > 1) {
    indexType numnodes = order - 1;
    numnodes *= numnodes;
    shape.resize(numnodes);
    shapeDerivative.resize(2, numnodes);
    indexType counter = 0;
    prec xV;
    prec eV;
    prec xVd;
    prec eVd;
    for (auto i = 2; i < order + 1; ++i) {
      LobattoShapes::getShape(xV, xVd, xsi, i);
      // LobattoShape(xV, xVd, xsi, i);
      for (auto j = 2; j < order + 1; ++j) {
        LobattoShapes::getShape(eV, eVd, eta, j);
        // LobattoShape(eV, eVd, eta, j);
        shape(counter) = eV * xV;
        shapeDerivative(0, counter) = xVd * eV;
        shapeDerivative(1, counter) = xV * eVd;
        ++counter;
      }
    }

  } else {
    shape.resize(0);
    shapeDerivative.resize(2, 0);
  }
}

auto LinearQuadrilateralData::getH1Shapes(indexType order,
                                          IntegrationPoint &IntegrationPt)
    -> H1Shapes {

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
      auto &tempEdge = *m_edges_pointers[0];
      prec orientation =
          tempEdge.getEdgeOrientation(this->m_verts[0], this->m_verts[1]);

      IntegrationPoint eInt;
      eInt.xi = xsi;
      auto eShape = tempEdge.getH1ShapesInternal(order, eInt);
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
      auto &tempEdge = *m_edges_pointers[1];
      prec orientation =
          tempEdge.getEdgeOrientation(this->m_verts[1], this->m_verts[2]);

      IntegrationPoint eInt;
      eInt.xi = eta;
      auto eShape = tempEdge.getH1ShapesInternal(order, eInt);
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
      auto &tempEdge = *m_edges_pointers[2];
      prec orientation =
          tempEdge.getEdgeOrientation(this->m_verts[3], this->m_verts[2]);

      IntegrationPoint eInt;
      eInt.xi = xsi;
      auto eShape = tempEdge.getH1ShapesInternal(order, eInt);
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
      auto &tempEdge = *m_edges_pointers[3];
      prec orientation =
          tempEdge.getEdgeOrientation(this->m_verts[0], this->m_verts[3]);

      IntegrationPoint eInt;
      eInt.xi = eta;
      auto eShape = tempEdge.getH1ShapesInternal(order, eInt);
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

auto LinearQuadrilateralData::getH1ShapesInternal(
    indexType order,
    IntegrationPoint &IntegrationPt, faceorientation orientation) -> H1Shapes {

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
void LinearQuadrilateralData::setHDivShapes(indexType meshId, indexType order,
                                            NodeTypes type) {
  indexType numNodes = order + 1;
  for (auto i = 0; i < 4; ++i) {
    auto &edgeTemp = *m_edges_pointers[i]; 
    edgeTemp.setNodeSet(meshId, numNodes, type);
  }

  if (order >= 1) {
    numNodes = 2 * order * order + 2 * order;
    this->setNodeSet(meshId, numNodes, type);
  }
}

void LinearQuadrilateralData::getHDivDofs(std::vector<DegreeOfFreedom *> &Dofs,
                                          indexType meshID, indexType order,
                                          NodeTypes type) {
  for (auto i = 0; i < 4; ++i) {
    auto &tempEdge = *m_edges_pointers[i];
    auto nodeList = tempEdge.getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }

  if (order >= 1) {
    auto nodeList = this->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }
}

void LinearQuadrilateralData::getHDivShapes(indexType order,
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
    auto &tempEdge = *m_edges_pointers[0];
    prec orientation =
        tempEdge.getEdgeOrientation(this->m_verts[0], this->m_verts[1]);
    prec orientationx =
        tempEdge.getEdgeOrientation(this->m_verts[0], this->m_verts[1]);
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
    auto &tempEdge = *m_edges_pointers[1];
    prec orientation =
        tempEdge.getEdgeOrientation(this->m_verts[1], this->m_verts[2]);
    prec orientationx =
        tempEdge.getEdgeOrientation(this->m_verts[1], this->m_verts[2]);
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
    auto &tempEdge = *m_edges_pointers[2];
    prec orientation =
        tempEdge.getEdgeOrientation(this->m_verts[2], this->m_verts[3]);
    prec orientationx =
        tempEdge.getEdgeOrientation(this->m_verts[3], this->m_verts[2]);
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
    auto &tempEdge = *m_edges_pointers[3];
    prec orientation =
        tempEdge.getEdgeOrientation(this->m_verts[3], this->m_verts[0]);
    prec orientationx =
        tempEdge.getEdgeOrientation(this->m_verts[0], this->m_verts[3]);
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

auto LinearQuadrilateralData::getHDivShapes(indexType order,
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
    auto &tempEdge = *m_edges_pointers[0];
    prec orientation =
        tempEdge.getEdgeOrientation(this->m_verts[0], this->m_verts[1]);
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
    auto &tempEdge = *m_edges_pointers[1];
    prec orientation =
        tempEdge.getEdgeOrientation(this->m_verts[1], this->m_verts[2]);
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
    auto &tempEdge = *m_edges_pointers[2];
    prec orientation =
        tempEdge.getEdgeOrientation(this->m_verts[2], this->m_verts[3]);
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
    auto &tempEdge = *m_edges_pointers[3];
    prec orientation =
        tempEdge.getEdgeOrientation(this->m_verts[3], this->m_verts[0]);
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

void LinearQuadrilateralData::setAllNodeBoundaryConditionMeshId(
    indexType meshId, indexType dof) {
  GeometryBaseData::setAllNodeBoundaryConditionMeshId(meshId, dof);

  for (auto &i : m_edges_pointers) {
    i->setAllNodeBoundaryConditionMeshId(meshId, dof);
  }
}


void LinearQuadrilateralData::setLoad(
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
        dx += shape.shapeDeriv(0, j) * m_verts_pointers[j]->getCoordinates();
        dy += shape.shapeDeriv(1, j) * m_verts_pointers[j]->getCoordinates();
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

void LinearQuadrilateralData::setPrescribedSolution(
    LoadList &loadlist, indexType meshid, ShapeFunctionTypes shapeType,
    indexType shapeOrder, Types::VectorX<prec> &Solution, indexType propNumber,
    Types::VectorX<prec> &direction, bool local, bool add) {
  if (shapeType == ShapeFunctionTypes::H1) {
    for (auto &i : m_verts_pointers) {
      i->setPrescribedSolution(loadlist, meshid, shapeType, shapeOrder,
                               Solution, propNumber, direction, local, add);
    }
  }
}

void LinearQuadrilateralData::setBoundaryCondition(
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

void LinearQuadrilateralData::setSpecialPlateShapes(indexType meshid,
                                                    indexType order,
                                                    NodeTypes type) {}
auto LinearQuadrilateralData::getSpecialPlateDofs(indexType meshID,
                                                  indexType order,
                                                  NodeTypes type)
    -> std::vector<DegreeOfFreedom *> {
  std::vector<DegreeOfFreedom *> Dofs;
  // auto &edge = pointers.getGeometryData()->getEdge(this->edges[0]);
  return Dofs;
}
auto LinearQuadrilateralData::getSpecialPlateShapes(IntegrationPoint &intPoint,
                                                    indexType order)
    -> SpecialPlateShapes {
  SpecialPlateShapes shapes;
  prec fValue;
  prec fDeriv;
  LobattoShapes::getShape(fValue, fDeriv, intPoint.xi, 0);
  return shapes;
}

void LinearQuadrilateralData::checkUpdateElement(EquationHandler &eqHandler,
                                                 GeometryData &geoData) {

  if (this->m_verts[0] == -1) {
    throw std::runtime_error("Element has no vertices");
  }

  const static std::array<indexType, 4> eeverts = {1, 2, 3, 0};

  for (auto i = 0; i < 4; ++i) {
    if (this->m_edges[i] == -1) {
      auto &Vert1 = geoData.getVertexData(this->m_verts[i]);

      auto edgeNums = Vert1.getConnectedEdges();
      bool edgeExists = false;
      bool search = true;
      if (edgeNums.size() == 0) {
        search = false;
      }
      indexType pos = 0;
      while (search) {
        auto &edge = geoData.getEdgeData(edgeNums[pos]);
        if (edge.hasVertices(this->m_verts[i], this->m_verts[eeverts[i]])) {
          edgeExists = true;
          search = false;
          this->m_edges[i] = edgeNums[pos];
        } else {
          pos++;
          if (pos >= static_cast<indexType>(edgeNums.size())) {
            search = false;
          }
        }
      }
      if (!edgeExists) {
        auto edgeNum = geoData.requestNewGeometryObject(
            eqHandler, GeometryTypes::LinearEdge);
        auto &edge = geoData.getEdgeData(edgeNum);
        std::vector<indexType> edgeVerts = {this->m_verts[i],
                                            this->m_verts[eeverts[i]]};
        edge.setVerts(geoData, edgeVerts);
        this->m_edges[i] = edgeNum;
      }
    }
  }
}

void LinearQuadrilateralData::flip() {

  std::reverse(this->m_verts.begin() + 1, this->m_verts.end());
  std::reverse(this->m_edges.begin(), this->m_edges.end());
  std::reverse(this->m_verts_pointers.begin() + 1,
               this->m_verts_pointers.end());
  std::reverse(this->m_edges_pointers.begin(), this->m_edges_pointers.end());
}

void LinearQuadrilateralData::rotate(indexType n) {
  while (n >= 4) {
    n -= 4;
  }
  while (n < 0) {
    n += 4;
  }
  std::rotate(this->m_verts.begin(), this->m_verts.begin() + n,
              this->m_verts.end());
  std::rotate(this->m_edges.begin(), this->m_edges.begin() + n,
              this->m_edges.end());

  std::rotate(this->m_verts_pointers.begin(),
              this->m_verts_pointers.begin() + n, this->m_verts_pointers.end());
  std::rotate(this->m_edges_pointers.begin(),
              this->m_edges_pointers.begin() + n, this->m_edges_pointers.end());
}

auto LinearQuadrilateralData::computeMeanCoordinate()
    -> Types::Vector3<prec> {

  Types::Vector3<prec> meanCoor = Types::Vector3<prec>::Zero();

  for (auto &i : this->m_verts_pointers) {
    meanCoor += i->getCoordinates();
  }
  meanCoor *= prec(0.25);
  return meanCoor;
}

auto LinearQuadrilateralData::getVertexNumber(indexType localNumber)
    -> indexType {
  return m_verts[localNumber];
}



const GeometryTypes LinearQuadrilateralData::type =
    GeometryTypes::LinearQuadrilateral;
} // namespace HierAMuS::Geometry
