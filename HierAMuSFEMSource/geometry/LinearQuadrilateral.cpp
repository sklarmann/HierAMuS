// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "MatrixTypes.h"
#include "equations/DegreeOfFreedom.h"
#include "geometry/Base.h"
#include "geometry/Faces.h"
#include "geometry/GeometryTypes.h"
#include "plot/vtkplotClassBase.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "shapefunctions/LegendreShapes.h"
#include "shapefunctions/LobattoShapes.h"
#include <exception>


#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <geometry/LinearQuadrilateral.h>
#include <pointercollection/pointercollection.h>
#include <sstream>

#include <equations/NodeSet.h>
#include <geometry/Edges.h>
#include <geometry/GeometryData.h>
#include <geometry/Vertex.h>
#include "geometry/GeometryData.h"

#include <loads/LoadList.h>

#include <stdexcept>
#include <vector>
#include <vtkCellType.h>

#include "HelperFunctions.h"

namespace HierAMuS::Geometry {

LinearQuadrilateral::LinearQuadrilateral()
    : verts({-1, -1, -1, -1}), edges({-1, -1, -1, -1}){};

LinearQuadrilateral::~LinearQuadrilateral() = default;

auto LinearQuadrilateral::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::LinearQuadrilateral::type;
}

void LinearQuadrilateral::setVerts(GeometryData &geoData, std::vector<indexType> &vertsIn) {
  if (vertsIn.size() == 4) {
    for (auto i = 0; i < 4; ++i) {
      this->verts[i] = vertsIn[i];
      auto &V1 = geoData.getVertex(vertsIn[i]);
      V1.connectFace(this->id);
    }
  } else {
    std::stringstream temp;
    temp << "Linear Quadrilateral geometry element has exactly 4 vertices. "
            "Cannot "
         << "add the given amount of " << vertsIn.size() << " vertices!";
    std::runtime_error(temp.str());
  }
}

void LinearQuadrilateral::setEdges(const std::vector<indexType> &edgesIn) {
  if (edgesIn.size() == 4) {
    for (auto i = 0; i < 4; ++i) {
      this->edges[i] = edgesIn[i];
    }
  } else {
    std::stringstream temp;
    temp << "Quadratic / Linear Quadrilateral geometry element has exactly 4 Edges. Cannot "
         << "add the given amount of " << edgesIn.size() << " edges!";
    std::runtime_error(temp.str());
  }
}

void LinearQuadrilateral::getVerts(std::vector<indexType> &vertsOut) {
  vertsOut.resize(4);
  for (auto i = 0; i < 4; ++i) {
    vertsOut[i] = this->verts[i];
  }
}

void LinearQuadrilateral::getVerts(PointerCollection &pointers,
                                   std::vector<Base *> &vertsOut) {
  vertsOut.clear();
  for (auto i = 0; i < 4; i++) {
    vertsOut.push_back(&pointers.getGeometryData()->getVertex(this->verts[i]));
  }
}

auto LinearQuadrilateral::getVertex(PointerCollection &pointers, indexType local_number)
    -> Geometry::Vertex * {
  if (local_number > 3) {
    throw std::runtime_error(
        "Requested a non-existing vertex from a linear / quadratic quadrilateral element.");
  }
  return &pointers.getGeometryData()->getVertex(this->verts[local_number]);
}

auto LinearQuadrilateral::getEdge(PointerCollection &pointers, indexType local_number) -> Geometry::Edges * {
  if (local_number > 3) {
    throw std::runtime_error(
        "Requested a non-existing edge from a linear quadrilateral element.");
  }
  return &pointers.getGeometryData()->getEdge(this->edges[local_number]);
}

void LinearQuadrilateral::getEdges(std::vector<indexType> &edgesOut) {
  edgesOut.resize(4);
  // edgesOut.insert(edgesOut.begin(), this->edges, this->edges + 4);
  for (auto i = 0; i < 4; ++i) {
    edgesOut[i] = this->edges[i];
  }
}

void LinearQuadrilateral::getEdges(PointerCollection &pointers, std::vector<Base *> &edgesOut) {
  edgesOut.clear();
  for (auto i = 0; i < 4; i++) {
    edgesOut.push_back(&pointers.getGeometryData()->getEdge(this->edges[i]));
  }
}

auto LinearQuadrilateral::getEdge(indexType local_number) -> indexType {
  return edges[local_number];
}

auto LinearQuadrilateral::hasVertices(indexType v1, indexType v2, indexType v3) -> bool {
  if(contains(this->verts, v1)){
    if(contains(this->verts, v2)){
      if(contains(this->verts, v3)){
        return true;
      }
    }
  }
  return false;
}

inline void LinearQuadrilateral::print(PointerCollection &pointers) {
  auto &Logger = pointers.getSPDLogger();

  Logger.debug("Linear Quadrilateral Face id: {:d}", this->id);
  Logger.debug("Vertices: {}", fmt::join(this->verts, " "));
  Logger.debug("Edges: {}", fmt::join(this->edges, " "));

  this->printEqInfo(pointers);
}

auto LinearQuadrilateral::getCoordinates(PointerCollection &pointers, prec xi,
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
    t1 = pointers.getGeometryData()->getVertex(this->verts[i]).getCoordinates();
    ret += t1 * vals[i];
  }
  return ret;
}

auto LinearQuadrilateral::getCoordinates(PointerCollection &pointers,
                                         IntegrationPoint &integrationPoint)
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
    t1 = pointers.getGeometryData()->getVertex(this->verts[i]).getCoordinates();
    ret += t1 * vals[i];
  }
  return ret;
}

auto LinearQuadrilateral::getTangent_G1(PointerCollection &pointers,
                                        IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  //auto shapes = this->getH1Shapes(pointers, 1, integrationPoint);
  Types::Vector3<prec> ret;
  //ret.setZero();
  //for (auto i = 0; i < 4; ++i) {
  //  ret +=
  //      shapes.shapeDeriv(0, i) *
  //      pointers.getGeometryData()->getVertex(this->verts[i]).getCoordinates();
  //}

  auto geoData = pointers.getGeometryData();
  auto [x1, dx1] = LobattoShapes::getShape(integrationPoint.xi, 0);
  auto [y1, dy1] = LobattoShapes::getShape(integrationPoint.eta, 0);
  auto [x2, dx2] = LobattoShapes::getShape(integrationPoint.xi, 1);
  auto [y2, dy2] = LobattoShapes::getShape(integrationPoint.eta, 1);

  ret = geoData->getVertex(this->verts[0]).getCoordinates() * (dx1 * y1);
  ret += geoData->getVertex(this->verts[1]).getCoordinates() * (dx2 * y1);
  ret += geoData->getVertex(this->verts[2]).getCoordinates() * (dx2 * y2);
  ret += geoData->getVertex(this->verts[3]).getCoordinates() * (dx1 * y2);

  return ret;
}

auto LinearQuadrilateral::getTangent_G2(PointerCollection &pointers,
                                        IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
    //auto shapes = this->getH1Shapes(pointers, 1, integrationPoint);
    Types::Vector3<prec> ret;
    //ret.setZero();
    //for (auto i = 0; i < 4; ++i) {
    //  ret +=
    //      shapes.shapeDeriv(1, i) *
    //      pointers.getGeometryData()->getVertex(this->verts[i]).getCoordinates();
    //}

  
  auto geoData = pointers.getGeometryData();
  auto [x1, dx1] = LobattoShapes::getShape(integrationPoint.xi, 0);
  auto [y1, dy1] = LobattoShapes::getShape(integrationPoint.eta, 0);
  auto [x2, dx2] = LobattoShapes::getShape(integrationPoint.xi, 1);
  auto [y2, dy2] = LobattoShapes::getShape(integrationPoint.eta, 1);

  ret = geoData->getVertex(this->verts[0]).getCoordinates() *  (x1 * dy1);
  ret += geoData->getVertex(this->verts[1]).getCoordinates() * (x2 * dy1);
  ret += geoData->getVertex(this->verts[2]).getCoordinates() * (x2 * dy2);
  ret += geoData->getVertex(this->verts[3]).getCoordinates() * (x1 * dy2);
  return ret;
}

auto LinearQuadrilateral::getFaceNormal(PointerCollection &pointers)
    -> Types::Vector3<prec> {
  auto &V1 = pointers.getGeometryData()->getVertex(this->verts[0]);
  auto &V2 = pointers.getGeometryData()->getVertex(this->verts[1]);
  auto &V3 = pointers.getGeometryData()->getVertex(this->verts[3]);

  Types::Vector3<prec> dx = V2.getCoordinates() - V1.getCoordinates();
  Types::Vector3<prec> dy = V3.getCoordinates() - V1.getCoordinates();

  Types::Vector3<prec> n = dx.cross(dy).normalized();
  return n;
}

auto LinearQuadrilateral::getOrientation(PointerCollection &pointers,
                                         indexType vertex1, indexType vertex2)
    -> faceorientation {
  faceorientation result = faceorientation::p_1;
  // look for first vertex position
  indexType i = 0;
  if (this->verts[0] == vertex1) {
    i = 0;
  } else if (this->verts[1] == vertex1) {
    i = 1;
  } else if (this->verts[2] == vertex1) {
    i = 2;
  } else if (this->verts[3] == vertex1) {
    i = 3;
  } else {
    throw std::runtime_error(
        "Vertex1 not found in faceorientation of linearquadrilateral");
  }

  // look for second vertex position
  indexType nextVertex = i + 1;
  indexType prevVertex = i - 1;
  nextVertex > 3 ? nextVertex = 0 : nextVertex = nextVertex;
  prevVertex < 0 ? prevVertex = 3 : prevVertex = prevVertex;
  if (this->verts[nextVertex] == vertex2) {
    if (i == 0) {
      result = faceorientation::p_1;
    } else if (i == 1) {
      result = faceorientation::p_2;
    } else if (i == 2) {
      result = faceorientation::p_3;
    } else if (i == 3) {
      result = faceorientation::p_4;
    }
  } else if (this->verts[prevVertex] == vertex2) {
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

void LinearQuadrilateral::modifyIntegrationpoint(IntegrationPoint &IP,
                                                 prec &shapeFactor,
                                                 faceorientation orientation) {
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

auto LinearQuadrilateral::getIntegrationPoints(PointerCollection &pointers, indexType elementId)
-> IntegrationPoints {
  IntegrationPoints temp =
      HierAMuS::PointerCollection::getIntegrationPoints(elementId);
  temp.setType(IntegrationType::Gauss2D);
  return temp;
}

auto LinearQuadrilateral::getJacobian(PointerCollection &pointers,
                                      IntegrationPoint &IntegrationPt)
    -> Types::MatrixXX<prec> {
  Types::MatrixXX<prec> jacobi;
  jacobi.resize(2, 2);
  //jacobi.setZero();

  prec N1x = IntegrationPt.eta * prec(0.25) - prec(0.25);
  prec N2x = prec(0.25) - IntegrationPt.eta * prec(0.25);
  prec N3x = IntegrationPt.eta * prec(0.25) + prec(0.25);
  prec N4x = -IntegrationPt.eta * prec(0.25) - prec(0.25);

  prec N1e = IntegrationPt.xi * prec(0.25) - prec(0.25);
  prec N2e = -IntegrationPt.xi * prec(0.25) - prec(0.25);
  prec N3e = IntegrationPt.xi * prec(0.25) + prec(0.25);
  prec N4e = prec(0.25) - IntegrationPt.xi * prec(0.25);

  // Vertex 1
  auto coord =
      pointers.getGeometryData()->getVertex(this->verts[0]).getCoordinates();
  jacobi(0, 0) = N1x * coord(0);
  jacobi(1, 0) = N1x * coord(1);
  jacobi(0, 1) = N1e * coord(0);
  jacobi(1, 1) = N1e * coord(1);
  
  // Vertex 2
  coord =
      pointers.getGeometryData()->getVertex(this->verts[1]).getCoordinates();
  jacobi(0, 0) += N2x * coord(0);
  jacobi(1, 0) += N2x * coord(1);
  jacobi(0, 1) += N2e * coord(0);
  jacobi(1, 1) += N2e * coord(1);
  
  // Vertex 3
  coord =
      pointers.getGeometryData()->getVertex(this->verts[2]).getCoordinates();
  jacobi(0, 0) += N3x * coord(0);
  jacobi(1, 0) += N3x * coord(1);
  jacobi(0, 1) += N3e * coord(0);
  jacobi(1, 1) += N3e * coord(1);
  
  // Vertex 4
  coord =
      pointers.getGeometryData()->getVertex(this->verts[3]).getCoordinates();
  jacobi(0, 0) += N4x * coord(0);
  jacobi(1, 0) += N4x * coord(1);
  jacobi(0, 1) += N4e * coord(0);
  jacobi(1, 1) += N4e * coord(1);
  
  return jacobi;
}

// H1Shapes

void LinearQuadrilateral::setH1Shapes(PointerCollection &pointers,
                                      indexType meshId, indexType order,
                                      NodeTypes type) {
  for (auto i = 0; i < 4; ++i) {
    auto &edgeTemp = pointers.getGeometryData()->getEdge(this->edges[i]);
    edgeTemp.setH1Shapes(pointers, meshId, order, type);
  }
  this->setH1ShapesInternal(pointers, meshId, order, type);
}

void LinearQuadrilateral::setH1ShapesInternal(PointerCollection &pointers,
                                              indexType meshId, indexType order,
                                              NodeTypes type) {
  if (order > 1) {
    indexType numNodes = order - 1;
    numNodes *= numNodes;

    this->setNodeSet(pointers, meshId, numNodes, type);
  }
}

auto LinearQuadrilateral::getH1Dofs(PointerCollection &pointers,
                                    indexType meshID, indexType order)
    -> std::vector<DegreeOfFreedom *> {
  std::vector<DegreeOfFreedom *> Dofs;
  this->getH1Dofs(pointers, Dofs, meshID, order);
  return Dofs;
}

void LinearQuadrilateral::getH1Dofs(PointerCollection &pointers,
                                    std::vector<DegreeOfFreedom *> &Dofs,
                                    indexType meshID, indexType order) {
  std::vector<DegreeOfFreedom *> tdofs;
  NodeSet *tempSet;
  for (auto i = 0; i < 4; ++i) {
    auto &tempVert = pointers.getGeometryData()->getVertex(this->verts[i]);
    tempSet = tempVert.getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
  if (order > 1) {
    for (auto i = 0; i < 4; ++i) {
      auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[i]);
      tempEdge.getH1DofsInternal(pointers, Dofs, meshID, order);
    }
    this->getH1DofsInternal(pointers, Dofs, meshID, order);
  }
}

void LinearQuadrilateral::getH1DofsInternal(
    PointerCollection &pointers, std::vector<DegreeOfFreedom *> &Dofs,
    indexType meshID, indexType order) {
  if (order > 1) {
    std::vector<DegreeOfFreedom *> tdofs;
    NodeSet *tempSet;
    tempSet = this->getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
}

auto LinearQuadrilateral::getH1Nodes(PointerCollection &pointers,
                                     indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  for (auto i : this->verts) {
    auto &V = pointers.getGeometryData()->getVertex(i);
    auto tempNodes = V.getNodesOfSet(pointers, meshID);
    nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
  }
  if (order > 1) {
    for (auto i : this->edges) {
      auto &E = pointers.getGeometryData()->getEdge(i);
      auto tempNodes = E.getH1NodesInternal(pointers, meshID, order);
      nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
    }
    auto tempNodes = this->getH1NodesInternal(pointers, meshID, order);
    nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
  }
  return nodes;
}

auto LinearQuadrilateral::getH1NodesInternal(PointerCollection &pointers,
                                             indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  if (order > 1) {
    indexType totnodes = (order - 1) * (order - 1);
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

auto LinearQuadrilateral::getHDivNodes(PointerCollection &pointers,
                                       indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  return {};
}

auto LinearQuadrilateral::getHDivNodesInternal(PointerCollection &pointers,
                                               indexType meshID,
                                               indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  if (order > 0) {
    indexType totnodes = (order - 1) * (order - 1);
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
  return {};
}

void LinearQuadrilateral::getH1Shapes(PointerCollection &pointers,
                                      indexType order,
                                      Types::VectorX<prec> &shape,
                                      Types::Matrix2X<prec> &shapeDerivative,
                                      prec xsi, prec eta) {

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
      auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[0]);
      prec orientation =
          tempEdge.getEdgeOrientation(this->verts[0], this->verts[1]);
      prec txsi = xsi * orientation;
      tempEdge.getH1ShapesInternal(pointers, order, tempshape,
                                   tempEdgeDerivative, txsi);
      tempEdgeDerivative *= orientation;
      for (auto i = 0; i < tempshape.rows(); ++i) {
        shape(counter) = tempshape(i) * s3;
        shapeDerivative(0, counter) = tempEdgeDerivative(i) * s3;
        shapeDerivative(1, counter) = tempshape(i) * ds3;
        ++counter;
      }
    }
    // Edge 2
    {
      auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[1]);
      prec orientation =
          tempEdge.getEdgeOrientation(this->verts[1], this->verts[2]);
      prec txsi = eta * orientation;
      tempEdge.getH1ShapesInternal(pointers, order, tempshape,
                                   tempEdgeDerivative, txsi);
      tempEdgeDerivative *= orientation;
      for (auto i = 0; i < tempshape.rows(); ++i) {
        shape(counter) = tempshape(i) * s2;
        shapeDerivative(1, counter) = tempEdgeDerivative(i) * s2;
        shapeDerivative(0, counter) = tempshape(i) * ds2;
        ++counter;
      }
    }
    // Edge 3
    {
      auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[2]);
      prec orientation =
          tempEdge.getEdgeOrientation(this->verts[3], this->verts[2]);
      prec txsi = xsi * orientation;
      tempEdge.getH1ShapesInternal(pointers, order, tempshape,
                                   tempEdgeDerivative, txsi);
      tempEdgeDerivative *= orientation;
      for (auto i = 0; i < tempshape.rows(); ++i) {
        shape(counter) = tempshape(i) * s4;
        shapeDerivative(0, counter) = tempEdgeDerivative(i) * s4;
        shapeDerivative(1, counter) = tempshape(i) * ds4;
        ++counter;
      }
    }
    // Edge 4
    {
      auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[3]);
      prec orientation =
          tempEdge.getEdgeOrientation(this->verts[0], this->verts[3]);
      prec txsi = eta * orientation;
      tempEdge.getH1ShapesInternal(pointers, order, tempshape,
                                   tempEdgeDerivative, txsi);
      tempEdgeDerivative *= orientation;
      for (auto i = 0; i < tempshape.rows(); ++i) {
        shape(counter) = tempshape(i) * s1;
        shapeDerivative(1, counter) = tempEdgeDerivative(i) * s1;
        shapeDerivative(0, counter) = tempshape(i) * ds1;
        ++counter;
      }
    }
    this->getH1ShapesInternal(pointers, order, tempshape, tempshapeDerivative,
                              xsi, eta);
    for (auto i = 0; i < tempshape.rows(); ++i) {
      shape(counter) = tempshape(i);
      shapeDerivative(0, counter) = tempshapeDerivative(0, i);
      shapeDerivative(1, counter) = tempshapeDerivative(1, i);
      ++counter;
    }
  }
}

void LinearQuadrilateral::getH1ShapesInternal(
    PointerCollection &pointers, indexType order, Types::VectorX<prec> &shape,
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

auto LinearQuadrilateral::getH1Shapes(PointerCollection &pointers,
                                      indexType order,
                                      IntegrationPoint &IntegrationPt)
    -> H1Shapes {

  H1Shapes shapes;
  indexType vertshapes = 4;
  indexType edgeShapes = (order - 1) * 4;
  indexType faceShapes = (order - 1) * (order - 1);

  indexType numshapes = vertshapes + edgeShapes + faceShapes;
  shapes.shapes.resize(numshapes);
  shapes.shapeDeriv.resize(2, numshapes);

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
      auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[0]);
      prec orientation =
          tempEdge.getEdgeOrientation(this->verts[0], this->verts[1]);

      IntegrationPoint eInt;
      eInt.xi = xsi * orientation;
      auto eShape = tempEdge.getH1ShapesInternal(pointers, order, eInt);
      eShape.shapeDeriv *= orientation;

     
      for (auto i = 0; i < eShape.shapes.rows(); ++i) {
        shapes.shapes(counter) = eShape.shapes(i) * s3;
        shapes.shapeDeriv(0, counter) = eShape.shapeDeriv(0,i) * s3;
        shapes.shapeDeriv(1, counter) = eShape.shapes(i) * ds3;
        ++counter;
      }
    }
    // Edge 2
    {
      auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[1]);
      prec orientation =
          tempEdge.getEdgeOrientation(this->verts[1], this->verts[2]);

      IntegrationPoint eInt;
      eInt.xi = eta * orientation;
      auto eShape = tempEdge.getH1ShapesInternal(pointers, order, eInt);
      eShape.shapeDeriv *= orientation;

      for (auto i = 0; i < eShape.shapes.rows(); ++i) {
        shapes.shapes(counter) = eShape.shapes(i) * s2;
        shapes.shapeDeriv(1, counter) = eShape.shapeDeriv(0,i) * s2;
        shapes.shapeDeriv(0, counter) = eShape.shapes(i) * ds2;
        ++counter;
      }
    }
    // Edge 3
    {
      auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[2]);
      prec orientation =
          tempEdge.getEdgeOrientation(this->verts[3], this->verts[2]);

      IntegrationPoint eInt;
      eInt.xi = xsi * orientation;
      auto eShape = tempEdge.getH1ShapesInternal(pointers, order, eInt);
      eShape.shapeDeriv *= orientation;
      
      for (auto i = 0; i < eShape.shapes.rows(); ++i) {
        shapes.shapes(counter) = eShape.shapes(i) * s4;
        shapes.shapeDeriv(0, counter) = eShape.shapeDeriv(0,i) * s4;
        shapes.shapeDeriv(1, counter) = eShape.shapes(i) * ds4;
        ++counter;
      }
    }
    // Edge 4
    {
      auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[3]);
      prec orientation =
          tempEdge.getEdgeOrientation(this->verts[0], this->verts[3]);

      IntegrationPoint eInt;
      eInt.xi = eta * orientation;
      auto eShape = tempEdge.getH1ShapesInternal(pointers, order, eInt);
      eShape.shapeDeriv *= orientation;
      
      for (auto i = 0; i < eShape.shapes.rows(); ++i) {
        shapes.shapes(counter) = eShape.shapes(i) * s1;
        shapes.shapeDeriv(1, counter) = eShape.shapeDeriv(0,i) * s1;
        shapes.shapeDeriv(0, counter) = eShape.shapes(i) * ds1;
        ++counter;
      }
    }

    auto inShapes = this->getH1ShapesInternal(pointers, order, IntegrationPt);
    for (auto i = 0; i < inShapes.shapes.rows(); ++i) {
      shapes.shapes(counter) = inShapes.shapes(i);
      shapes.shapeDeriv(0, counter) = inShapes.shapeDeriv(0, i);
      shapes.shapeDeriv(1, counter) = inShapes.shapeDeriv(1, i);
      ++counter;
    }
  }

  return shapes;
}

auto LinearQuadrilateral::getH1ShapesInternal(PointerCollection &pointers,
                                              indexType order,
                                              IntegrationPoint &IntegrationPt,
                                              faceorientation orientation)
    -> H1Shapes {
  H1Shapes shapes;

  prec xsi = IntegrationPt.xi;
  prec eta = IntegrationPt.eta;

  indexType numShapes = order - 1;
  numShapes *= numShapes;
  shapes.shapes.resize(numShapes);
  shapes.shapeDeriv.resize(2, numShapes);
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
void LinearQuadrilateral::setHDivShapes(PointerCollection &pointers,
                                        indexType meshId,
                                        indexType order,
                                        NodeTypes type) {
  indexType numNodes = order + 1;
  for (auto i = 0; i < 4; ++i) {
    auto &edgeTemp = pointers.getGeometryData()->getEdge(this->edges[i]);
    edgeTemp.setNodeSet(pointers, meshId, numNodes, type);
  }

  if (order >= 1) {
    numNodes = 2 * order * order + 2 * order;
    this->setNodeSet(pointers, meshId, numNodes, type);
  }
}

void LinearQuadrilateral::getHDivDofs(PointerCollection &pointers,
                                      std::vector<DegreeOfFreedom *> &Dofs,
                                      indexType meshID, indexType order,
                                      NodeTypes type) {
  std::vector<DegreeOfFreedom *> tdofs;
  NodeSet *tempSet;
  for (auto i = 0; i < 4; ++i) {
    auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[i]);
    tempSet = tempEdge.getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }

  if (order >= 1) {
    tempSet = this->getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
}

void LinearQuadrilateral::getHDivShapes(PointerCollection &pointers,
                                        indexType order,
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
    auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[0]);
    prec orientation =
        tempEdge.getEdgeOrientation(this->verts[0], this->verts[1]);
    prec orientationx =
        tempEdge.getEdgeOrientation(this->verts[0], this->verts[1]);
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
    auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[1]);
    prec orientation =
        tempEdge.getEdgeOrientation(this->verts[1], this->verts[2]);
    prec orientationx =
        tempEdge.getEdgeOrientation(this->verts[1], this->verts[2]);
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
    auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[2]);
    prec orientation =
        tempEdge.getEdgeOrientation(this->verts[2], this->verts[3]);
    prec orientationx =
        tempEdge.getEdgeOrientation(this->verts[3], this->verts[2]);
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
    auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[3]);
    prec orientation =
        tempEdge.getEdgeOrientation(this->verts[3], this->verts[0]);
    prec orientationx =
        tempEdge.getEdgeOrientation(this->verts[0], this->verts[3]);
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

auto LinearQuadrilateral::getHDivShapes(PointerCollection &pointers,
                                        indexType order,
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
    auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[0]);
    prec orientation =
        tempEdge.getEdgeOrientation(this->verts[0], this->verts[1]);
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
    auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[1]);
    prec orientation =
        tempEdge.getEdgeOrientation(this->verts[1], this->verts[2]);
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
    auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[2]);
    prec orientation =
        tempEdge.getEdgeOrientation(this->verts[2], this->verts[3]);
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
    auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[3]);
    prec orientation =
        tempEdge.getEdgeOrientation(this->verts[3], this->verts[0]);
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

void LinearQuadrilateral::setAllNodeBoundaryConditionMeshId(
    PointerCollection &pointers, indexType meshId, indexType dof) {
  Base::setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);

  for (auto i = 0; i < 4; i++) {
    auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[i]);
    tempEdge.setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);
  }
}

void LinearQuadrilateral::geometryToParaview(PointerCollection &pointers,
                                             vtkPlotInterface &paraviewAdapter,
                                             indexType mainMesh,
                                             indexType subMesh) {
  indexType numPoints = 4;
  std::vector<indexType> points(numPoints);
  points.clear();
  for (auto i : this->verts) {
    auto &vert = pointers.getGeometryData()->getVertex(i);
    vert.geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
    points.push_back(i);
  }

  paraviewAdapter.addCell(mainMesh, subMesh, this->id, 1, points, 4, VTK_QUAD);
}

void LinearQuadrilateral::computeWeightsParaview(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh) {

  auto GP = this->getIntegrationPoints(pointers,-1);
  GP.setOrder(2);

  for (auto i : GP) {
    auto jaco = this->getJacobian(pointers, i);
    auto shapes = this->getH1Shapes(pointers, 1, i);
    prec dA = jaco.determinant() * i.weight;

    for (auto j = 0; j < 4; ++j) {
      auto &vert = pointers.getGeometryData()->getVertex(this->verts[j]);
      std::vector<prec> val;
      val.push_back(shapes.shapes(j) * dA);
      paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, val, vert.getId(),
                                           1, paraviewNames::weightName());
    }
  }
}

void LinearQuadrilateral::H1SolutionToParaview(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType meshId, indexType order,
    std::string &name) {
  std::vector<DegreeOfFreedom *> Dofs;
  this->getH1Dofs(pointers, Dofs, meshId, order);
  auto solution = pointers.getSolutionState()->getSolution(Dofs);
  for (auto i = 0; i < 4; ++i) {
    auto &V = pointers.getGeometryData()->getVertex(this->verts[i]);
    std::vector<prec> sol(3);
    for (auto j = 0; j < 3; ++j) {
      sol[j] = solution(3 * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, V.getId(), sol, 3, name);
  }
}

void LinearQuadrilateral::H1DataToParaview(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, Types::VectorX<prec> &Data,
    indexType numberComponents, indexType order, std::string &name) {
  for (auto i = 0; i < 4; ++i) {
    auto &V = pointers.getGeometryData()->getVertex(this->verts[i]);
    std::vector<prec> sol(numberComponents);
    for (auto j = 0; j < numberComponents; ++j) {
      sol[j] = Data(numberComponents * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, V.getId(), sol,
                                 numberComponents, name);
  }
}

void LinearQuadrilateral::projectDataToParaviewVertices(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {

  auto shapes = this->getH1Shapes(pointers, 1, IntegrationPt);
  std::vector<prec> vals(numberComponents);

  auto jaco = this->getJacobian(pointers, IntegrationPt);
  auto dA = jaco.determinant() * IntegrationPt.weight;
  for (auto i = 0; i < 4; ++i) {
    auto &V = pointers.getGeometryData()->getVertex(this->verts[i]);
    for (auto j = 0; j < numberComponents; ++j) {
      vals[j] = data(j) * shapes.shapes(i) * dA;
    }
    paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, vals, V.getId(),
                                         numberComponents, name);
  }
}

void LinearQuadrilateral::setLoad(
    PointerCollection &pointers, indexType meshid, ShapeFunctionTypes shapeType,
    indexType shapeOrder, Types::VectorX<prec> &Loads, indexType propNumber,
    Types::VectorX<prec> &direction, bool local, bool add) {

  if (shapeType == ShapeFunctionTypes::H1) {
    auto GP = this->getIntegrationPoints(pointers,-1);
    GP.setOrder(shapeOrder * 2);
    auto nodes = this->getH1Nodes(pointers, meshid, shapeOrder);
    Types::VectorX<prec> tempLoads(nodes.size() * 3);
    tempLoads.setZero();
    for (auto i : GP) {
      auto shape = this->getH1Shapes(pointers, shapeOrder, i);
      Types::Vector3<prec> dx;
      Types::Vector3<prec> dy;
      dx.setZero();
      dy.setZero();
      for (auto j = 0; j < 4; ++j) {
        auto &vert = pointers.getGeometryData()->getVertex(this->verts[j]);
        dx += shape.shapeDeriv(0, j) * vert.getCoordinates();
        dy += shape.shapeDeriv(1, j) * vert.getCoordinates();
      }
      auto dA = dx.cross(dy).norm() * i.weight;
      indexType counter = 0;
      for (auto j = 0; j < nodes.size(); ++j) {
        auto tempDofs = nodes[j]->getDegreesOfFreedom(pointers);
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
      auto tempDofs = i->getDegreesOfFreedom(pointers);
      for (auto j : tempDofs) {
        pointers.getLoadList()->setLoad(propNumber, j->getId(),
                                        tempLoads(counter), add);
        counter++;
      }
    }
  }
}

void LinearQuadrilateral::setPrescribedSolution(
    PointerCollection &pointers, indexType meshid, ShapeFunctionTypes shapeType,
    indexType shapeOrder, Types::VectorX<prec> &Solution, indexType propNumber,
    Types::VectorX<prec> &direction, bool local, bool add)
{
  if (shapeType == ShapeFunctionTypes::H1)
  {
    for (auto i:this->verts)
    {
      auto &Vert = pointers.getGeometryData()->getVertex(i);
      Vert.setPrescribedSolution(
          pointers, meshid,
          shapeType, shapeOrder,
          Solution, propNumber,
          direction, local, add);
    }
  }
}

void LinearQuadrilateral::setBoundaryCondition(
    PointerCollection &pointers, indexType meshId, indexType order,
    ShapeFunctionTypes shapeType, Types::Vector3<indexType> &dofs, bool set) {
  if (shapeType == ShapeFunctionTypes::H1) {
    auto nodes = this->getH1Nodes(pointers, meshId, order);
    for (auto &i : nodes) {
      if (set) {
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            i->setBoundaryCondition(pointers, j);
          } else {
            i->unsetBoundaryCondition(pointers, j);
          }
        }
      } else {
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            i->setBoundaryCondition(pointers, j);
          }
        }
      }
    }
  }
}

void LinearQuadrilateral::setSpecialPlateShapes(PointerCollection &pointers,
                                                indexType meshid,
                                                indexType order,
                                                NodeTypes type) {}
auto LinearQuadrilateral::getSpecialPlateDofs(PointerCollection &pointers,
                                              indexType meshID, indexType order,
                                              NodeTypes type)
    -> std::vector<DegreeOfFreedom *> {
  std::vector<DegreeOfFreedom *> Dofs;
  // auto &edge = pointers.getGeometryData()->getEdge(this->edges[0]);
  return Dofs;
}
auto LinearQuadrilateral::getSpecialPlateShapes(PointerCollection &pointers,
                                                IntegrationPoint &intPoint,
                                                indexType order)
    -> SpecialPlateShapes {
  SpecialPlateShapes shapes;
  prec fValue;
  prec fDeriv;
  LobattoShapes::getShape(fValue, fDeriv, intPoint.xi, 0);
  return shapes;
}

void LinearQuadrilateral::checkUpdateElement(GeometryData &geoData) {

  if (this->verts[0] == -1) {
    throw std::runtime_error("Element has no vertices");
  }

  const static std::array<indexType,4> eeverts ={1,2,3,0};


  for (auto i=0;i<4;++i){
    if (this->edges[i] == -1) {
      auto &Vert1 = geoData.getVertex(this->verts[i]);

      auto edgeNums = Vert1.getConnectedEdges();
      bool edgeExists = false;
      bool search = true;
      if(edgeNums.size() == 0){
        search = false;
      }
      indexType pos = 0;
      while(search) {
        auto &edge = geoData.getEdge(edgeNums[pos]);
        if(edge.hasVertices(this->verts[i], this->verts[eeverts[i]])) {
          edgeExists = true;
          search = false;
          this->edges[i] = edgeNums[pos];
        }  else {
          pos++;
          if (pos >= edgeNums.size()) {
            search = false;
          }
        }
      }
      if (!edgeExists) {
        auto edgeNum = geoData.requestNewGeometryObject(GeometryTypes::LinearEdge);
        auto &edge = geoData.getEdge(edgeNum);
        std::vector<indexType> edgeVerts = {this->verts[i], this->verts[eeverts[i]]};
        edge.setVerts(geoData, edgeVerts);
        this->edges[i] = edgeNum;
      }
    }
  }
}

void LinearQuadrilateral::flip()
{

  std::reverse(this->verts.begin() + 1, this->verts.end());
  std::reverse(this->edges.begin(), this->edges.end());
}

void LinearQuadrilateral::rotate(indexType n)
{
  while (n>=4)
  {
    n -= 4;
  }
  while (n<0)
  {
    n += 4;
  }
  std::rotate(this->verts.begin(), this->verts.begin() + n, this->verts.end());
  std::rotate(this->edges.begin(), this->edges.begin() + n, this->edges.end());
  
}

auto LinearQuadrilateral::computeMeanCoordinate(PointerCollection &pointers)
    -> Types::Vector3<prec> {

  Types::Vector3<prec> meanCoor;
  meanCoor.setZero();

  for (auto i:this->verts)
  {
    meanCoor += pointers.getGeometryData()->getVertex(i).getCoordinates();
  }
  meanCoor *= prec(0.25);
  return meanCoor;
}

const GeometryTypes LinearQuadrilateral::type =
    GeometryTypes::LinearQuadrilateral;
} // namespace HierAMuS::Geometry
