// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "geometry/GeometryData.h"
#include "datatypes.h"
#include "geometry/GeometryBaseData.h"
#include <geometry/Faces/LinearTriangleData.h>

#include "plot/vtkplotClass.h"

#include <shapefunctions/KernelShapes.h>
#include <shapefunctions/LobattoShapes.h>

#include <exception>

#include "GenericNodes.h"
#include <geometry/Edges/EdgesData.h>
#include <geometry/VertexData.h>

#include "LoadList.h"

#include <iomanip>

#include <vtkCellType.h>

#include "HelperFunctions.h"

#include "shapefunctions/LegendreShapes.h"

#include "geometry/Faces/LinearTriangleRuntime.h"

#include "shapefunctions/IntegrationsPoints/helperClasses/IntegrationPointsManagement.h"

namespace HierAMuS::Geometry {

LinearTriangleData::LinearTriangleData() : FacesDataInterface(){};

LinearTriangleData::~LinearTriangleData() = default;

auto LinearTriangleData::getRuntimeObject(GeometryData &geoData)
    -> std::shared_ptr<FacesRuntime> {
  return std::make_shared<LinearTriangleRuntime>(geoData, *this);
}

auto LinearTriangleData::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::LinearTriangleData::type;
}

auto LinearTriangleData::getCoordinates(prec xi,
                                        prec eta) -> Types::Vector3<prec> {
  Types::Vector3<prec> ret;
  Types::Vector3<prec> t1;
  ret.setZero();
  prec s1;
  prec ds1;
  prec s2;
  prec ds2;
  prec s3;
  prec ds3;
  std::array<prec, 3> vals;

  LobattoShapes::getShape(s1, ds1, xi, 0);
  LobattoShapes::getShape(s2, ds2, xi, 1);
  LobattoShapes::getShape(s3, ds3, eta, 1);

  vals[0] = s1 - s3;
  vals[1] = s2;
  vals[2] = s3;

  for (auto i = 0; i < 3; ++i) {
    t1 = m_verts_pointers[i]->getCoordinates();
    ret += t1 * vals[i];
  }
  return ret;
}

auto LinearTriangleData::getCoordinates(IntegrationPoint &integrationPoint)
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
  std::array<prec, 3> vals;

  LobattoShapes::getShape(s1, ds1, integrationPoint.xi, 0);
  LobattoShapes::getShape(s2, ds2, integrationPoint.xi, 1);
  LobattoShapes::getShape(s3, ds3, integrationPoint.eta, 1);

  vals[0] = s1 - s3;
  vals[1] = s2;
  vals[2] = s3;

  for (auto i = 0; i < 3; ++i) {
    t1 = m_verts_pointers[i]->getCoordinates();
    ret += t1 * vals[i];
  }
  return ret;
}

void LinearTriangleData::checkUpdateElement(EquationHandler &eqHandler,
                                            GeometryData &geoData) {

  if (this->m_verts[0] == -1) {
    throw std::runtime_error("Element has no vertices");
  }

  const static std::array<indexType, 3> eeverts = {1, 2, 0};

  for (auto i = 0; i < 3; ++i) {
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

auto LinearTriangleData::getTangent_G1(IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  auto shapes = this->getH1Shapes(1, integrationPoint);
  Types::Vector3<prec> ret;
  ret.setZero();
  for (auto i = 0; i < 3; ++i) {
    ret += shapes.shapeDeriv(0, i) * m_verts_pointers[i]->getCoordinates();
  }
  return ret;
}

auto LinearTriangleData::getTangent_G2(IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  auto shapes = this->getH1Shapes(1, integrationPoint);
  Types::Vector3<prec> ret;
  ret.setZero();
  for (auto i = 0; i < 3; ++i) {
    ret += shapes.shapeDeriv(1, i) * m_verts_pointers[i]->getCoordinates();
  }
  return ret;
}

auto LinearTriangleData::getFaceNormal()
    -> Types::Vector3<prec> {

  auto &V1 = *m_verts_pointers[0];
  auto &V2 = *m_verts_pointers[1];
  auto &V3 = *m_verts_pointers[2];

  Types::Vector3<prec> dx = V2.getCoordinates() - V1.getCoordinates();
  Types::Vector3<prec> dy = V3.getCoordinates() - V1.getCoordinates();

  auto n = dx.cross(dy);
  n.normalize();
  return n;
}

auto LinearTriangleData::getOrientation(indexType vertex1, indexType vertex2)
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
  } else {
    throw std::runtime_error(
        "Vertex1 not found in faceorientation of LinearTriangle");
  }

  // look for second vertex position
  indexType nextVertex = i + 1;
  indexType prevVertex = i - 1;
  if (nextVertex > 2)
    nextVertex = 0;
  if (prevVertex < 0)
    prevVertex = 2;
  if (this->m_verts[nextVertex] == vertex2) {
    if (i == 0) {
      result = faceorientation::p_1;
    } else if (i == 1) {
      result = faceorientation::p_2;
    } else if (i == 2) {
      result = faceorientation::p_3;
    }
  } else if (this->m_verts[prevVertex] == vertex2) {
    if (i == 1) {
      result = faceorientation::n_1;
    } else if (i == 2) {
      result = faceorientation::n_2;
    } else if (i == 0) {
      result = faceorientation::n_4;
    }
  } else {
    throw std::runtime_error(
        "Something wrong in faceorientation for lineartriangle");
  }
  // result = faceorientation::p_1;
  return result;
}

auto LinearTriangleData::getIntegrationPoints(indexType elementId)
    -> IntegrationPoints {
  IntegrationPoints temp =
      IntegrationPointsManagement::getIntegrationsPoints(elementId);
  temp.setType(IntegrationType::Gauss2DTriangle);
  return temp;
}

auto LinearTriangleData::getJacobian(IntegrationPoint &IntegrationPt)
    -> Types::Matrix22<prec> {
  Types::Matrix22<prec> jacobi;
  jacobi.setZero();

  auto shapes = this->getH1Shapes(1, IntegrationPt);

  for (auto i = 0; i < 3; ++i) {
    auto coord = m_verts_pointers[i]->getCoordinates();
    jacobi(0, 0) += shapes.shapeDeriv(0, i) * coord(0);
    jacobi(0, 1) += shapes.shapeDeriv(1, i) * coord(0);
    jacobi(1, 0) += shapes.shapeDeriv(0, i) * coord(1);
    jacobi(1, 1) += shapes.shapeDeriv(1, i) * coord(1);
  }
  return jacobi;
}

// ----------------H1------------------

void LinearTriangleData::setH1Shapes(indexType meshId, indexType order,
                                     NodeTypes type) {
  for (auto i = 0; i < 3; ++i) {
    m_edges_pointers[i]->setH1Shapes(meshId, order, type);
  }
  this->setH1ShapesInternal(meshId, order, type);
}

void LinearTriangleData::setH1ShapesInternal(indexType meshId, indexType order,
                                             NodeTypes type) {
  if (order > 2) {
    indexType numNodes = (order - 1) * (order - 2) / 2;
    this->setNodeSet(meshId, numNodes, type);
  }
}

auto LinearTriangleData::getH1Dofs(indexType meshID, indexType order)
    -> std::vector<DegreeOfFreedom *> {
  std::vector<DegreeOfFreedom *> Dofs;
  this->getH1Dofs(Dofs, meshID, order);
  return Dofs;
}

void LinearTriangleData::getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs,
                                   indexType meshID, indexType order) {
  for (auto i = 0; i < 3; ++i) {
    auto nodeList = m_verts_pointers[i]->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }
  if (order > 1) {
    for (auto i = 0; i < 3; ++i) {
      m_edges_pointers[i]->getH1DofsInternal(Dofs, meshID, order);
    }
    this->getH1DofsInternal(Dofs, meshID, order);
  }
}

void LinearTriangleData::getH1DofsInternal(std::vector<DegreeOfFreedom *> &Dofs,
                                           indexType meshID, indexType order) {
  if (order > 2) {
    auto nodeList = this->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }
}

auto LinearTriangleData::getH1Nodes(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  for (auto i : this->m_verts) {
    auto tempNodes = m_verts_pointers[i]->getNodesOfSet(meshID);
    nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
  }
  if (order > 1) {
    for (auto i : this->m_edges) {
      auto tempNodes = m_edges_pointers[i]->getH1NodesInternal(meshID, order);
      nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
    }
    auto tempNodes = this->getH1NodesInternal(meshID, order);
    nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
  }
  return nodes;
}

auto LinearTriangleData::getH1NodesInternal(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  if (order > 2) {
    indexType totnodes = (order - 1) * (order - 2) / 2;
    auto tempnodes = this->getNodesOfSet(meshID);
    nodes.insert(nodes.end(), tempnodes.begin(), tempnodes.end());
    if (totnodes != tempnodes.size()) {
      std::stringstream ss;
      ss << "Error in LinearTriangle::getH1NodesInternal: "
         << "Expected " << totnodes << " nodes, got " << tempnodes.size()
         << " instead.";
      throw std::runtime_error(ss.str());
    }
  }
  return nodes;
}

auto LinearTriangleData::getH1NodesList(indexType meshID, indexType order)
    -> MeshIdNodeList {
  return MeshIdNodeList(meshID);
}

auto LinearTriangleData::getH1Shapes(indexType order,
                                     IntegrationPoint &IntegrationPt)
    -> H1Shapes {

  indexType vertshapes = 3;
  indexType edgeShapes = (order - 1) * 3;
  indexType faceShapes = 0;
  if (order > 2) {
    faceShapes = (order - 1) * (order - 2) / 2;
  }

  indexType numshapes = vertshapes + edgeShapes + faceShapes;

  H1Shapes shapes(numshapes, 2);

  prec xsi = IntegrationPt.xi;
  prec eta = IntegrationPt.eta;

  prec Lam1 = (prec(1) + eta) * prec(0.5);
  prec Lam2 = -(xsi + eta) * prec(0.5);
  prec Lam3 = (prec(1) + xsi) * prec(0.5);

  prec dLam2x = -prec(0.5);
  prec dLam3x = prec(0.5);

  prec dLam1e = prec(0.5);
  prec dLam2e = -prec(0.5);

  shapes.shapes(0) = Lam2;
  shapes.shapes(1) = Lam3;
  shapes.shapes(2) = Lam1;

  shapes.shapeDeriv(0, 0) = dLam2x;
  shapes.shapeDeriv(0, 1) = dLam3x;
  shapes.shapeDeriv(0, 2) = prec(0);

  shapes.shapeDeriv(1, 0) = dLam2e;
  shapes.shapeDeriv(1, 1) = prec(0);
  shapes.shapeDeriv(1, 2) = dLam1e;

  indexType counter = 3;
  if (order > 1) {
    prec tempshape;
    prec tempEdgeDerivative;

    // Edge 1
    {
      prec orientation = m_edges_pointers[0]->getEdgeOrientation(
          this->m_verts[0], this->m_verts[1]);
      // std::cout << "Face: " << this->id <<  "Edge: " << tempEdge.getId()
      //           << " orientation: " << orientation << std::endl;

      for (auto i = 0; i < order - 1; ++i) {
        KernelShapes::getShape(tempshape, tempEdgeDerivative, (Lam3 - Lam2), i);

        shapes.shapes(counter) = tempshape * Lam2 * Lam3;
        shapes.shapeDeriv(0, counter) =
            (dLam2x * Lam3 + Lam2 * dLam3x) * tempshape +
            Lam2 * Lam3 * tempEdgeDerivative * (dLam3x - dLam2x);
        shapes.shapeDeriv(1, counter) =
            dLam2e * Lam3 * tempshape +
            Lam2 * Lam3 * tempEdgeDerivative * (-dLam2e);
        if (i & 1) {
          shapes.shapes(counter) *= orientation;
          shapes.shapeDeriv(0, counter) *= orientation;
          shapes.shapeDeriv(1, counter) *= orientation;
        }
        ++counter;
      }
    }

    // Edge 2
    {
      prec orientation = m_edges_pointers[1]->getEdgeOrientation(
          this->m_verts[1], this->m_verts[2]);

      // std::cout << "Face: " << this->id << "Edge: " << tempEdge.getId()
      //           << " orientation: " << orientation << std::endl;

      for (auto i = 0; i < order - 1; ++i) {
        KernelShapes::getShape(tempshape, tempEdgeDerivative, (Lam1 - Lam3), i);

        shapes.shapes(counter) = tempshape * Lam3 * Lam1;
        shapes.shapeDeriv(0, counter) =
            Lam1 * dLam3x * tempshape +
            Lam1 * Lam3 * tempEdgeDerivative * (-dLam3x);
        shapes.shapeDeriv(1, counter) =
            dLam1e * Lam3 * tempshape +
            Lam1 * Lam3 * tempEdgeDerivative * (dLam1e);
        if (i & 1) {
          shapes.shapes(counter) *= orientation;
          shapes.shapeDeriv(0, counter) *= orientation;
          shapes.shapeDeriv(1, counter) *= orientation;
        }
        ++counter;
      }
    }
    // Edge 3
    {
      prec orientation = m_edges_pointers[2]->getEdgeOrientation(
          this->m_verts[2], this->m_verts[0]);
      // std::cout << "Face: " << this->id << "Edge: " << tempEdge.getId()
      //           << " orientation: " << orientation << std::endl;

      for (auto i = 0; i < order - 1; ++i) {
        KernelShapes::getShape(tempshape, tempEdgeDerivative, (Lam2 - Lam1), i);

        shapes.shapes(counter) = tempshape * Lam1 * Lam2;
        shapes.shapeDeriv(0, counter) =
            Lam1 * dLam2x * tempshape +
            Lam1 * Lam2 * tempEdgeDerivative * (dLam2x);
        shapes.shapeDeriv(1, counter) =
            (dLam1e * Lam2 + Lam1 * dLam2e) * tempshape +
            Lam1 * Lam2 * tempEdgeDerivative * (dLam2e - dLam1e);
        if (i & 1) {
          shapes.shapes(counter) *= orientation;
          shapes.shapeDeriv(0, counter) *= orientation;
          shapes.shapeDeriv(1, counter) *= orientation;
        }
        ++counter;
      }
    }
  }

  // Face
  if (order > 2) {
    prec tempshape1, tempshape2;
    prec tempShapeDerivative1, tempShapeDerivative2;
    Types::Vector2<prec> Lam_Deriv;

    Lam_Deriv(0) = Lam1 * dLam2x * Lam3 + Lam1 * Lam2 * dLam3x;
    Lam_Deriv(1) = dLam1e * Lam2 * Lam3 + Lam1 * dLam2e * Lam3;
    prec lams = Lam1 * Lam2 * Lam3;

    for (auto i = 0; i < order - 2; ++i) {
      for (auto j = 0; j < order - 2 - i; ++j) {
        KernelShapes::getShape(tempshape1, tempShapeDerivative1, Lam3 - Lam2,
                               i);
        KernelShapes::getShape(tempshape2, tempShapeDerivative2, Lam2 - Lam1,
                               j);

        shapes.shapes(counter) = lams * tempshape1 * tempshape2;
        shapes.shapeDeriv(0, counter) =
            Lam_Deriv(0) * tempshape1 * tempshape2 +
            lams * tempShapeDerivative1 * (dLam3x - dLam2x) * tempshape2 +
            lams * tempshape1 * tempShapeDerivative2 * (dLam2x);
        shapes.shapeDeriv(1, counter) =
            Lam_Deriv(1) * tempshape1 * tempshape2 +
            lams * tempShapeDerivative1 * (-dLam2e) * tempshape2 +
            lams * tempshape1 * tempShapeDerivative2 * (dLam2e - dLam1e);
        ++counter;
      }
    }
  }

  return shapes;
}


// ----------------HDiv------------------

auto LinearTriangleData::getHDivNodes(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  return {};
}

auto LinearTriangleData::getHDivNodesInternal(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  if (order > 1) {
    indexType totnodes = (order - 1) * (order - 1);
    auto tempnodes = this->getNodesOfSet(meshID);
    nodes.insert(nodes.end(), tempnodes.begin(), tempnodes.end());
    if (totnodes != tempnodes.size()) {
      std::stringstream ss;
      ss << "Error in LinearTriangle::getH1NodesInternal: "
         << "Expected " << totnodes << " nodes, got " << tempnodes.size()
         << " instead.";
      throw std::runtime_error(ss.str());
    }
  }
  return nodes;
  return {};
}

void LinearTriangleData::setHDivShapes(indexType meshId, indexType order,
                                       NodeTypes type) {
  indexType numNodes = order + 1;
  for (auto i = 0; i < 3; ++i) {
    m_edges_pointers[i]->setNodeSet(meshId, numNodes, type);
  }
  if (order > 1) {
    numNodes = order * order - 1;
    this->setNodeSet(meshId, numNodes, type);
  }
}

void LinearTriangleData::getHDivDofs(std::vector<DegreeOfFreedom *> &Dofs,
                                     indexType meshID, indexType order,
                                     NodeTypes type) {
  for (auto i = 0; i < 3; ++i) {
    auto nodeList = m_edges_pointers[i]->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }
  if (order > 1) {
    auto nodeList = this->getNodeSetNodeListMeshId(meshID);
    nodeList.addDofsToVector(Dofs);
  }
}

auto LinearTriangleData::getHDivShapes(indexType order,
                                       IntegrationPoint &IntegrationPt)
    -> HDivShapes {
  HDivShapes shapes;

  prec xsi = IntegrationPt.xi;
  prec eta = IntegrationPt.eta;

  indexType edgeShapes = (order + 1) * 3;
  indexType faceShapes = 0;
  if (order > 1) {
    faceShapes = order * order - 1;
  }

  indexType numshapes = edgeShapes + faceShapes;
  shapes.shapes.resize(2, numshapes);
  shapes.shapeDeriv.resize(1, numshapes);
  shapes.shapeDeriv.setZero();

  prec Lk1, dLk1, Lk2, dLk2;
  prec gam0x, gam0e, gam1x, gam1e;
  prec Lam1, Lam2, Lam3;

  Lam1 = (eta + prec(1.0)) / prec(2.0);
  Lam2 = (-xsi - eta) / prec(2);
  Lam3 = (xsi + prec(1.0)) / prec(2.0);

  int counter = 0;

  // edge1
  {
    prec orientation = m_edges_pointers[0]->getEdgeOrientation(
        this->m_verts[0], this->m_verts[1]);
    gam0x = Lam3 * orientation;
    gam0e = (-Lam3 - Lam2) * orientation;
    gam1x = Lam3;
    gam1e = (Lam2 - Lam3);

    shapes.shapes(0, counter) = gam0x;
    shapes.shapes(1, counter) = gam0e;
    ++counter;
    if (order > 0) {
      shapes.shapes(0, counter) = gam1x;
      shapes.shapes(1, counter) = gam1e;
      ++counter;
      if (order > 1) {
        for (auto i = 2; i <= order; i++) {
          LegendreShapes::getShape(Lk1, dLk1, orientation * (Lam3 - Lam2),
                                   i - 1);
          LegendreShapes::getShape(Lk2, dLk2, orientation * (Lam3 - Lam2),
                                   i - 2);

          shapes.shapes(0, counter) =
              (prec(2.0) * prec(i) - prec(1.0)) / prec(i) * Lk1 * gam1x -
              (prec(i) - prec(1.0)) / prec(i) * Lk2 * gam0x;
          shapes.shapes(1, counter) =
              (prec(2.0) * prec(i) - prec(1.0)) / prec(i) * Lk1 * gam1e -
              (prec(i) - prec(1.0)) / prec(i) * Lk2 * gam0e;
          ++counter;
        }
      }
    }
  }

  // edge2
  {
    prec orientation = m_edges_pointers[1]->getEdgeOrientation(
        this->m_verts[1], this->m_verts[2]);
    // gam0x = Lam3 * 1/sqrt(2);
    // gam0e = Lam1 * 1/sqrt(2);
    // gam1x = -Lam3 * 1/sqrt(2);
    // gam1e = Lam1 * 1/sqrt(2);
    gam0x = Lam3 * orientation;
    gam0e = Lam1 * orientation;
    gam1x = -Lam3;
    gam1e = Lam1;

    shapes.shapes(0, counter) = gam0x;
    shapes.shapes(1, counter) = gam0e;
    ++counter;
    if (order > 0) {
      shapes.shapes(0, counter) = gam1x;
      shapes.shapes(1, counter) = gam1e;
      ++counter;
      if (order > 1) {
        for (auto i = 2; i <= order; i++) {
          LegendreShapes::getShape(Lk1, dLk1, orientation * (Lam1 - Lam3),
                                   i - 1);
          LegendreShapes::getShape(Lk2, dLk2, orientation * (Lam1 - Lam3),
                                   i - 2);

          shapes.shapes(0, counter) =
              (prec(2.0) * prec(i) - prec(1.0)) / prec(i) * Lk1 * gam1x -
              (prec(i) - prec(1.0)) / prec(i) * Lk2 * gam0x;
          shapes.shapes(1, counter) =
              (prec(2.0) * prec(i) - prec(1.0)) / prec(i) * Lk1 * gam1e -
              (prec(i) - prec(1.0)) / prec(i) * Lk2 * gam0e;
          ++counter;
        }
      }
    }
  }

  // edge3
  {
    prec orientation = m_edges_pointers[2]->getEdgeOrientation(
        this->m_verts[2], this->m_verts[0]);
    gam0x = (-Lam2 - Lam1) * orientation;
    gam0e = Lam1 * orientation;
    gam1x = (-Lam2 + Lam1);
    gam1e = -Lam1;

    shapes.shapes(0, counter) = gam0x;
    shapes.shapes(1, counter) = gam0e;
    ++counter;
    if (order > 0) {
      shapes.shapes(0, counter) = gam1x;
      shapes.shapes(1, counter) = gam1e;
      ++counter;
      if (order > 1) {
        for (auto i = 2; i <= order; i++) {
          LegendreShapes::getShape(Lk1, dLk1, orientation * (Lam2 - Lam1),
                                   i - 1);
          LegendreShapes::getShape(Lk2, dLk2, orientation * (Lam2 - Lam1),
                                   i - 2);

          shapes.shapes(0, counter) =
              (prec(2.0) * prec(i) - prec(1.0)) / prec(i) * Lk1 * gam1x -
              (prec(i) - prec(1.0)) / prec(i) * Lk2 * gam0x;
          shapes.shapes(1, counter) =
              (prec(2.0) * prec(i) - prec(1.0)) / prec(i) * Lk1 * gam1e -
              (prec(i) - prec(1.0)) / prec(i) * Lk2 * gam0e;
          ++counter;
        }
      }
    }
  }

  // Face Edge-base
  if (order > 1) {
    for (auto i = 2; i <= order; i++) {
      LegendreShapes::getShape(Lk2, dLk2, Lam3 - Lam2, i - 2);
      shapes.shapes(0, counter) = Lam3 * Lam2 * Lk2;
      shapes.shapes(1, counter) = 0;
      ++counter;
      LegendreShapes::getShape(Lk2, dLk2, Lam1 - Lam3, i - 2);
      shapes.shapes(0, counter) = Lam1 * Lam3 * Lk2 * -prec(1.0);
      shapes.shapes(1, counter) = Lam1 * Lam3 * Lk2 * prec(1.0);
      ++counter;
      LegendreShapes::getShape(Lk2, dLk2, Lam2 - Lam1, i - 2);
      shapes.shapes(0, counter) = 0;
      shapes.shapes(1, counter) = Lam1 * Lam2 * Lk2 * -prec(1.0);
      ++counter;
    }
  }
  // Face genuine
  if (order > 2) {
    for (auto i = 0; i < order - 2; i++) {
      for (auto j = 0; j < order - 2 - i; j++) {
        LegendreShapes::getShape(Lk1, dLk1, Lam3 - Lam2, i);
        LegendreShapes::getShape(Lk2, dLk2, Lam2 - Lam1, j);
        shapes.shapes(0, counter) = Lam1 * Lam2 * Lam3 * Lk1 * Lk2;
        shapes.shapes(1, counter) = 0;
        ++counter;
        shapes.shapes(0, counter) = 0;
        shapes.shapes(1, counter) = Lam1 * Lam2 * Lam3 * Lk1 * Lk2;
        ++counter;
      }
    }
  }

  return shapes;
}

void LinearTriangleData::setAllNodeBoundaryConditionMeshId(indexType meshId,
                                                           indexType dof) {
  GeometryBaseData::setAllNodeBoundaryConditionMeshId(meshId, dof);

  for (auto &i : m_edges_pointers) {
    i->setAllNodeBoundaryConditionMeshId(meshId, dof);
  }
}

void LinearTriangleData::setLoad(
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
      for (auto j = 0; j < 3; ++j) {
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

void LinearTriangleData::setBoundaryCondition(indexType meshId, indexType order,
                                              ShapeFunctionTypes shapeType,
                                              Types::Vector3<indexType> &dofs,
                                              bool set) {
  if (shapeType == ShapeFunctionTypes::H1) {
    auto nodes = this->getH1Nodes(meshId, order);
    for (auto i : nodes) {
      if (set) {
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            i->setBoundaryCondition(j);
          } else {
            i->unsetBoundaryCondition(dofs(j));
          }
        }
      } else {
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            i->setBoundaryCondition(dofs(j));
          }
        }
      }
    }
  } else if (shapeType == ShapeFunctionTypes::HDiv) {
    for (auto i = 0; i < 3; ++i) {
      if (dofs[i] != 0) {
        GeometryBaseData::setAllNodeBoundaryConditionMeshId(meshId, i);
      }
    }
  }
}


void LinearTriangleData::flip() {
  throw std::runtime_error(
      "ERROR in LinearTriangle::flip: Method not implemented!");
}

void LinearTriangleData::rotate(indexType n) {
  throw std::runtime_error(
      "ERROR in LinearTriangle::rotate: Method not implemented!");
}

auto LinearTriangleData::computeMeanCoordinate()
    -> Types::Vector3<prec> {
  throw std::runtime_error("ERROR in LinearTriangle::computeMeanCoordinate: "
                           "Method not implemented!");
  return Types::Vector3<prec>();
}

const GeometryTypes LinearTriangleData::type = GeometryTypes::LinearTriangle;

} // namespace HierAMuS::Geometry
