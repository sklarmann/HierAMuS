// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "geometry/GeometryData.h"
#include "datatypes.h"
#include "geometry/GeometryBaseData.h"
#include <geometry/Faces/LinearTriangleRuntime.h>

#include "plot/vtkplotClass.h"

#include <shapefunctions/KernelShapes.h>
#include <shapefunctions/LobattoShapes.h>

#include <exception>

#include "geometry/Edges/EdgesRuntime.h"
#include "geometry/Edges/EdgeH1ShapesInterface.h"
#include "geometry/Faces/LinearTriangleData.h"
#include "geometry/VertexRuntime.h"
#include <geometry/Edges/EdgesData.h>
#include <geometry/VertexData.h>

#include "LoadList.h"

#include <iomanip>

#include <vtkCellType.h>

#include "HelperFunctions.h"

#include "shapefunctions/LegendreShapes.h"

#include "DegreeOfFreedom.h"
#include "GenericNodes.h"

#include "shapefunctions/IntegrationsPoints/helperClasses/IntegrationPointsManagement.h"

namespace HierAMuS::Geometry {

LinearTriangleRuntime::LinearTriangleRuntime(GeometryData &geoData,
                                             LinearTriangleData &base_element)
    : FacesRuntimeDataInterface(geoData, base_element){};

LinearTriangleRuntime::~LinearTriangleRuntime() = default;

auto LinearTriangleRuntime::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::LinearTriangleRuntime::type;
}

void LinearTriangleRuntime::print(spdlog::logger &Log) {
  m_Face_Data_Element.print(Log);
}

void LinearTriangleRuntime::checkUpdateElement(EquationHandler &eqHandler,
                                               GeometryData &geoData) {
  m_Face_Data_Element.checkUpdateElement(eqHandler, geoData);
}

auto LinearTriangleRuntime::getTangent_G1(IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  auto shapes = this->getH1Shapes(1, integrationPoint);
  Types::Vector3<prec> ret;
  ret.setZero();
  for (auto i = 0; i < 3; ++i) {
    ret += shapes.shapeDeriv(0, i) * m_Vertices[i]->getCoordinates();
  }
  return ret;
}

auto LinearTriangleRuntime::getTangent_G2(IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  auto shapes = this->getH1Shapes(1, integrationPoint);
  Types::Vector3<prec> ret;
  ret.setZero();
  for (auto i = 0; i < 3; ++i) {
    ret += shapes.shapeDeriv(1, i) * m_Vertices[i]->getCoordinates();
  }
  return ret;
}

auto LinearTriangleRuntime::getFaceNormal()
    -> Types::Vector3<prec> {

  auto &V1 = m_Vertices[0];
  auto &V2 = m_Vertices[1];
  auto &V3 = m_Vertices[2];

  Types::Vector3<prec> dx = V2->getCoordinates() - V1->getCoordinates();
  Types::Vector3<prec> dy = V3->getCoordinates() - V1->getCoordinates();

  auto n = dx.cross(dy);
  n.normalize();
  return n;
}

auto LinearTriangleRuntime::getOrientation(indexType vertex1, indexType vertex2)
    -> faceorientation {
  return m_Face_Data_Element.getOrientation(vertex1, vertex2);
}

auto LinearTriangleRuntime::getIntegrationPoints(indexType elementId)
    -> IntegrationPoints {
  IntegrationPoints temp =
      IntegrationPointsManagement::getIntegrationsPoints(elementId);
  temp.setType(IntegrationType::Gauss2DTriangle);
  return temp;
}

auto LinearTriangleRuntime::getJacobian(IntegrationPoint &IntegrationPt)
    -> Types::Matrix22<prec> {
  Types::Matrix22<prec> jacobi;
  jacobi.setZero();

  auto shapes = this->getH1Shapes(1, IntegrationPt);

  for (auto i = 0; i < 3; ++i) {
    auto coord = m_Vertices[i]->getCoordinates();
    jacobi(0, 0) += shapes.shapeDeriv(0, i) * coord(0);
    jacobi(0, 1) += shapes.shapeDeriv(1, i) * coord(0);
    jacobi(1, 0) += shapes.shapeDeriv(0, i) * coord(1);
    jacobi(1, 1) += shapes.shapeDeriv(1, i) * coord(1);
  }
  return jacobi;
}

// ----------------H1------------------

void LinearTriangleRuntime::setH1Shapes(indexType meshId, indexType order,
                                        NodeTypes type) {
  for (auto i = 0; i < 3; ++i) {
    m_Edges[i]->getH1Edge()->setH1Shapes(meshId, order, type);
  }
  this->setH1ShapesInternal(meshId, order, type);
}

void LinearTriangleRuntime::setH1ShapesInternal(indexType meshId,
                                                indexType order,
                                                NodeTypes type) {
  if (order > 2) {
    indexType numNodes = (order - 1) * (order - 2) / 2;
    this->setNodeSet(meshId, numNodes, type);
  }
}

auto LinearTriangleRuntime::getH1Dofs(indexType meshID, indexType order)
    -> std::vector<DegreeOfFreedom *> {
  std::vector<DegreeOfFreedom *> Dofs;
  this->getH1Dofs(Dofs, meshID, order);
  return Dofs;
}

void LinearTriangleRuntime::getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs,
                                      indexType meshID, indexType order) {
  for (auto i = 0; i < 3; ++i) {
    auto nodeList = m_Vertices[i]->getNodeSetNodeListMeshId(meshID);
    auto tdofs = nodeList.getDegreesOfFreedom();
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
  if (order > 1) {
    for (auto i = 0; i < 3; ++i) {
      m_Edges[i]->getH1Edge()->getH1DofsInternal(Dofs, meshID, order);
    }
    this->getH1DofsInternal(Dofs, meshID, order);
  }
}

void LinearTriangleRuntime::getH1DofsInternal(
    std::vector<DegreeOfFreedom *> &Dofs,
    indexType meshID, indexType order) {
  if (order > 2) {
    auto nodeList = this->getNodeSetNodeListMeshId(meshID);
    auto tdofs = nodeList.getDegreesOfFreedom();
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
}

auto LinearTriangleRuntime::getH1Nodes(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  for (auto &i : m_Vertices) {
    auto setNodes = i->getNodeSetNodeListMeshId(meshID);
    for (auto &j : setNodes) {
      nodes.push_back(&j);
    }
  }
  if (order > 1) {
    for (auto &i : m_Edges) {
      auto tempNodes = i->getH1Edge()->getH1NodesInternal(meshID, order);
      nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
    }
    auto tempNodes = this->getH1NodesInternal(meshID, order);
    nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
  }
  return nodes;
}

auto LinearTriangleRuntime::getH1NodesInternal(indexType meshID,
                                               indexType order)
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

auto LinearTriangleRuntime::getH1NodesList(indexType meshID, indexType order)
    -> MeshIdNodeList {
  return MeshIdNodeList(meshID);
}

auto LinearTriangleRuntime::getH1Shapes(indexType order,
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
      prec orientation = m_Edges[0]->getEdgeOrientation(m_Vertices[0]->getId(),
                                                        m_Vertices[1]->getId());
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
      prec orientation = m_Edges[1]->getEdgeOrientation(m_Vertices[1]->getId(),
                                                        m_Vertices[2]->getId());

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
      prec orientation = m_Edges[2]->getEdgeOrientation(m_Vertices[2]->getId(),
                                                        m_Vertices[0]->getId());
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

auto LinearTriangleRuntime::getH1ShapesInternal(indexType order,
                                                IntegrationPoint &IntegrationPt,
                                                faceorientation orientation)
    -> H1Shapes  {
  return H1Shapes(0,0);
}

// ----------------HDiv------------------

auto LinearTriangleRuntime::getHDivNodes(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  return {};
}

auto LinearTriangleRuntime::getHDivNodesInternal(indexType meshID,
                                                 indexType order)
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

void LinearTriangleRuntime::setHDivShapes(indexType meshId, indexType order,
                                          NodeTypes type) {
  indexType numNodes = order + 1;
  for (auto i = 0; i < 3; ++i) {
    m_Edges[i]->setNodeSet(meshId, numNodes, type);
  }
  if (order > 1) {
    numNodes = order * order - 1;
    this->setNodeSet(meshId, numNodes, type);
  }
}

void LinearTriangleRuntime::getHDivDofs(std::vector<DegreeOfFreedom *> &Dofs,
                                        indexType meshID, indexType order,
                                        NodeTypes type) {
  for (auto i = 0; i < 3; ++i) {
    auto nodeList = m_Edges[i]->getNodeSetNodeListMeshId(meshID);
    auto tdofs(nodeList.getDegreesOfFreedom());
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
  if (order > 1) {
    auto nodeList = this->getNodeSetNodeListMeshId(meshID);
    auto tdofs = nodeList.getDegreesOfFreedom();
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
}

auto LinearTriangleRuntime::getHDivShapes(indexType order,
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
    prec orientation = m_Edges[0]->getEdgeOrientation(m_Vertices[0]->getId(),
                                                      m_Vertices[1]->getId());
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
    prec orientation = m_Edges[1]->getEdgeOrientation(m_Vertices[1]->getId(),
                                                      m_Vertices[2]->getId());
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
    prec orientation = m_Edges[2]->getEdgeOrientation(m_Vertices[2]->getId(),
                                                      m_Vertices[0]->getId());
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

void LinearTriangleRuntime::getHDivShapes(indexType order,
                                          Types::Matrix2X<prec> &shape,
                                          Types::VectorX<prec> &dshape, prec xi,
                                          prec eta) {}

void LinearTriangleRuntime::setAllNodeBoundaryConditionMeshId(indexType meshId,
                                                              indexType dof) {
  m_Face_Data_Element.setAllNodeBoundaryConditionMeshId(meshId, dof);
}

void LinearTriangleRuntime::setLoad(
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

void LinearTriangleRuntime::setBoundaryCondition(
    indexType meshId, indexType order, ShapeFunctionTypes shapeType,
    Types::Vector3<indexType> &dofs, bool set) {
  m_Face_Data_Element.setBoundaryCondition(meshId, order, shapeType, dofs, set);
}

void LinearTriangleRuntime::geometryToParaview(
    vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh) {
  indexType numPoints = 3;
  std::vector<indexType> points(numPoints);
  points.clear();
  for (auto &i : m_Vertices) {
    i->geometryToParaview(paraviewAdapter, mainMesh, subMesh);
    points.push_back(i->getId());
  }

  paraviewAdapter.addCell(mainMesh, subMesh, this->m_Face_Data_Element.getId(),
                          1, points, 3, VTK_TRIANGLE);
}

void LinearTriangleRuntime::computeWeightsParaview(
    vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh) {

  auto GP = this->getIntegrationPoints(-1);
  GP.setOrder(2);

  for (auto i : GP) {
    auto jaco = this->getJacobian(i);
    auto shapes = this->getH1Shapes(1, i);
    prec dA = jaco.determinant() * i.weight;

    for (auto i = 0; i < 3; ++i) {
      auto &vert = m_Vertices[i];
      std::vector<prec> val;
      val.push_back(shapes.shapes(i) * dA);
      paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, val,
                                           vert->getId(), 1,
                                           paraviewNames::weightName());
    }
  }
}

void LinearTriangleRuntime::H1SolutionToParaview(
    vtkPlotInterface &paraviewAdapter, indexType mainMesh, indexType subMesh,
    indexType order, Types::VectorX<prec> &solution, std::string &name) {
  for (auto i = 0; i < 3; ++i) {
    std::vector<prec> sol(3);
    for (auto j = 0; j < 3; ++j) {
      sol[j] = solution(3 * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, m_Vertices[i]->getId(), sol,
                                 3, name);
  }
}

void LinearTriangleRuntime::H1DataToParaview(
    vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, Types::VectorX<prec> &Data,
    indexType numberComponents, indexType order, std::string &name) {
  for (auto i = 0; i < 3; ++i) {
    std::vector<prec> sol(numberComponents);
    for (auto j = 0; j < numberComponents; ++j) {
      sol[j] = Data(numberComponents * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, m_Vertices[i]->getId(), sol,
                                 numberComponents, name);
  }
}

void LinearTriangleRuntime::projectDataToParaviewVertices(
    vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {

  auto shapes = this->getH1Shapes(1, IntegrationPt);
  std::vector<prec> vals(numberComponents);

  auto jaco = this->getJacobian(IntegrationPt);
  auto dA = jaco.determinant() * IntegrationPt.weight;
  for (auto i = 0; i < 3; ++i) {
    for (auto j = 0; j < numberComponents; ++j) {
      vals[j] = data(j) * shapes.shapes(i) * dA;
    }
    paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, vals,
                                         m_Vertices[i]->getId(),
                                         numberComponents, name);
  }
}

void LinearTriangleRuntime::flip() {
  throw std::runtime_error(
      "ERROR in LinearTriangle::flip: Method not implemented!");
}

void LinearTriangleRuntime::rotate(indexType n) {
  throw std::runtime_error(
      "ERROR in LinearTriangle::rotate: Method not implemented!");
}

auto LinearTriangleRuntime::computeMeanCoordinate()
    -> Types::Vector3<prec> {
  throw std::runtime_error("ERROR in LinearTriangle::computeMeanCoordinate: "
                           "Method not implemented!");
  return Types::Vector3<prec>();
}

void LinearTriangleRuntime::set_geometry_pointers(GeometryData &geoData) {
  m_Face_Data_Element.set_geometry_pointers(geoData);
}

const GeometryTypes LinearTriangleRuntime::type = GeometryTypes::LinearTriangle;

} // namespace HierAMuS::Geometry
