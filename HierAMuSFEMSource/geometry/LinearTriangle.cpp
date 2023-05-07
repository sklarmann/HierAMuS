// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "geometry/Base.h"
#include <geometry/LinearTriangle.h>

#include "equations/DegreeOfFreedom.h"

#include "geometry/GeometryData.h"
#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <pointercollection/pointercollection.h>
#include <shapefunctions/KernelShapes.h>
#include <shapefunctions/LobattoShapes.h>

#include <exception>

#include "geometry/GeometryData.h"
#include <equations/NodeSet.h>
#include <geometry/Edges.h>
#include <geometry/GeometryData.h>
#include <geometry/Vertex.h>

#include <loads/LoadList.h>

#include <iomanip>

#include <vtkCellType.h>

#include "HelperFunctions.h"

namespace HierAMuS::Geometry {

LinearTriangle::LinearTriangle() : verts({-1, -1, -1}), edges({-1, -1, -1}){};

LinearTriangle::~LinearTriangle() = default;

auto LinearTriangle::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::LinearTriangle::type;
}

void LinearTriangle::setVerts(GeometryData &geoData,
                              std::vector<indexType> &vertsIn) {
  if (vertsIn.size() == 3) {
    for (auto i = 0; i < 3; ++i) {
      this->verts[i] = vertsIn[i];
      auto &V1 = geoData.getVertex(vertsIn[i]);
      V1.connectFace(this->id);
    }
  } else {
    std::stringstream temp;
    temp << "Linear Triangle geometry element has exactly 3 vertices. "
            "Cannot "
         << "add the given amount of " << vertsIn.size() << " vertices!";
    std::runtime_error(temp.str());
  }
}

inline void LinearTriangle::print(PointerCollection &pointers) {
  auto &Logger = pointers.getSPDLogger();

  Logger.debug("Linear Triangle Face id:  {:>12}", this->id);
  Logger.debug("Vertices:  {}", fmt::join(verts," "));
  Logger.debug("Edges:     {}", fmt::join(edges," "));

  this->printEqInfo(pointers);
}

void LinearTriangle::setEdges(const std::vector<indexType> &edgesIn) {
  if (edgesIn.size() == 3) {
    for (auto i = 0; i < 3; ++i) {
      this->edges[i] = edgesIn[i];
    }
  } else {
    std::stringstream temp;
    temp << "Linear Triangle geometry element has exactly 3 Edges. Cannot "
         << "add the given amount of " << edgesIn.size() << " edges!";
    std::runtime_error(temp.str());
  }
}

void LinearTriangle::getVerts(std::vector<indexType> &vertsOut) {
  vertsOut.clear();
  for (auto &i : this->verts) {
    vertsOut.push_back(i);
  }
}

void LinearTriangle::getVerts(PointerCollection &pointers,
                              std::vector<Base *> &vertsOut) {
  vertsOut.clear();
  for (auto &i : this->verts) {
    vertsOut.push_back(&pointers.getGeometryData()->getVertex(i));
  }
}

auto LinearTriangle::getVertex(PointerCollection &pointers,
                               indexType local_number) -> Geometry::Vertex * {
  if (local_number > 2) {
    throw std::runtime_error(
        "Requested a non-existing vertex from a linear triangle element.");
  }
  return &pointers.getGeometryData()->getVertex(this->verts[local_number]);
}

auto LinearTriangle::getEdge(PointerCollection &pointers,
                             indexType local_number) -> Geometry::Edges * {
  if (local_number > 2) {
    throw std::runtime_error(
        "Requested a non-existing edge from a linear triangle element.");
  }
  return &pointers.getGeometryData()->getEdge(this->edges[local_number]);
}

void LinearTriangle::getEdges(std::vector<indexType> &edgesOut) {
  edgesOut.clear();
  for (auto &i : this->edges) {
    edgesOut.push_back(i);
  }
}

auto LinearTriangle::getEdge(indexType local_number) -> indexType {
  return edges[local_number];
}

auto LinearTriangle::hasVertices(indexType v1, indexType v2, indexType v3)
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

auto LinearTriangle::getCoordinates(PointerCollection &pointers, prec xi,
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
    t1 = pointers.getGeometryData()->getVertex(this->verts[i]).getCoordinates();
    ret += t1 * vals[i];
  }
  return ret;
}

auto LinearTriangle::getCoordinates(PointerCollection &pointers,
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
  std::array<prec, 3> vals;

  LobattoShapes::getShape(s1, ds1, integrationPoint.xi, 0);
  LobattoShapes::getShape(s2, ds2, integrationPoint.xi, 1);
  LobattoShapes::getShape(s3, ds3, integrationPoint.eta, 1);

  vals[0] = s1 - s3;
  vals[1] = s2;
  vals[2] = s3;

  for (auto i = 0; i < 3; ++i) {
    t1 = pointers.getGeometryData()->getVertex(this->verts[i]).getCoordinates();
    ret += t1 * vals[i];
  }
  return ret;
}

void LinearTriangle::checkUpdateElement(GeometryData &geoData) {

  if (this->verts[0] == -1) {
    throw std::runtime_error("Element has no vertices");
  }

  const static std::array<indexType, 3> eeverts = {1, 2, 0};

  for (auto i = 0; i < 3; ++i) {
    if (this->edges[i] == -1) {
      auto &Vert1 = geoData.getVertex(this->verts[i]);

      auto edgeNums = Vert1.getConnectedEdges();
      bool edgeExists = false;
      bool search = true;
      if (edgeNums.size() == 0) {
        search = false;
      }
      indexType pos = 0;
      while (search) {
        auto &edge = geoData.getEdge(edgeNums[pos]);
        if (edge.hasVertices(this->verts[i], this->verts[eeverts[i]])) {
          edgeExists = true;
          search = false;
          this->edges[i] = edgeNums[pos];
        } else {
          pos++;
          if (pos >= edgeNums.size()) {
            search = false;
          }
        }
      }
      if (!edgeExists) {
        auto edgeNum =
            geoData.requestNewGeometryObject(GeometryTypes::LinearEdge);
        auto &edge = geoData.getEdge(edgeNum);
        std::vector<indexType> edgeVerts = {this->verts[i],
                                            this->verts[eeverts[i]]};
        edge.setVerts(geoData, edgeVerts);
        this->edges[i] = edgeNum;
      }
    }
  }
}

auto LinearTriangle::getTangent_G1(PointerCollection &pointers,
                                   IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  auto shapes = this->getH1Shapes(pointers, 1, integrationPoint);
  Types::Vector3<prec> ret;
  ret.setZero();
  for (auto i = 0; i < 3; ++i) {
    ret +=
        shapes.shapeDeriv(0, i) *
        pointers.getGeometryData()->getVertex(this->verts[i]).getCoordinates();
  }
  return ret;
}

auto LinearTriangle::getTangent_G2(PointerCollection &pointers,
                                   IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  auto shapes = this->getH1Shapes(pointers, 1, integrationPoint);
  Types::Vector3<prec> ret;
  ret.setZero();
  for (auto i = 0; i < 3; ++i) {
    ret +=
        shapes.shapeDeriv(1, i) *
        pointers.getGeometryData()->getVertex(this->verts[i]).getCoordinates();
  }
  return ret;
}

auto LinearTriangle::getFaceNormal(PointerCollection &pointers)
    -> Types::Vector3<prec> {
  auto &V1 = pointers.getGeometryData()->getVertex(this->verts[0]);
  auto &V2 = pointers.getGeometryData()->getVertex(this->verts[1]);
  auto &V3 = pointers.getGeometryData()->getVertex(this->verts[2]);

  Types::Vector3<prec> dx = V2.getCoordinates() - V1.getCoordinates();
  Types::Vector3<prec> dy = V3.getCoordinates() - V1.getCoordinates();

  auto n = dx.cross(dy);
  n.normalize();
  return n;
}

auto LinearTriangle::getOrientation(PointerCollection &pointers,
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
  } else {
    throw std::runtime_error(
        "Vertex1 not found in faceorientation of LinearTriangle");
  }

  // look for second vertex position
  indexType nextVertex = i + 1;
  indexType prevVertex = i - 1;
  nextVertex > 2 ? nextVertex = 0 : nextVertex = nextVertex;
  prevVertex < 0 ? prevVertex = 2 : prevVertex = prevVertex;
  if (this->verts[nextVertex] == vertex2) {
    if (i == 0) {
      result = faceorientation::p_1;
    } else if (i == 1) {
      result = faceorientation::p_2;
    } else if (i == 2) {
      result = faceorientation::p_3;
    }
  } else if (this->verts[prevVertex] == vertex2) {
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

auto LinearTriangle::getIntegrationPoints(PointerCollection &pointers, indexType elementId)
-> IntegrationPoints {
  IntegrationPoints temp =
      HierAMuS::PointerCollection::getIntegrationPoints(elementId);
  temp.setType(IntegrationType::Gauss2DTriangle);
  return temp;
}

auto LinearTriangle::getJacobian(PointerCollection &pointers,
                                 IntegrationPoint &IntegrationPt)
    -> Types::MatrixXX<prec> {
  Types::MatrixXX<prec> jacobi;
  jacobi.resize(2, 2);
  jacobi.setZero();

  auto shapes = this->getH1Shapes(pointers, 1, IntegrationPt);

  for (auto i = 0; i < 3; ++i) {
    auto coord =
        pointers.getGeometryData()->getVertex(this->verts[i]).getCoordinates();
    jacobi(0, 0) += shapes.shapeDeriv(0, i) * coord(0);
    jacobi(0, 1) += shapes.shapeDeriv(1, i) * coord(0);
    jacobi(1, 0) += shapes.shapeDeriv(0, i) * coord(1);
    jacobi(1, 1) += shapes.shapeDeriv(1, i) * coord(1);
  }
  return jacobi;
}

// ----------------H1------------------

void LinearTriangle::setH1Shapes(PointerCollection &pointers, indexType meshId,
                                 indexType order, NodeTypes type) {
  for (auto i = 0; i < 3; ++i) {
    auto &edgeTemp = pointers.getGeometryData()->getEdge(this->edges[i]);
    edgeTemp.setH1Shapes(pointers, meshId, order, type);
  }
  this->setH1ShapesInternal(pointers, meshId, order, type);
}

void LinearTriangle::setH1ShapesInternal(PointerCollection &pointers,
                                         indexType meshId, indexType order,
                                         NodeTypes type) {
  if (order > 2) {
    indexType numNodes = (order - 1) * (order - 2) / 2;
    this->setNodeSet(pointers, meshId, numNodes, type);
  }
}

auto LinearTriangle::getH1Dofs(PointerCollection &pointers, indexType meshID,
                               indexType order)
    -> std::vector<DegreeOfFreedom *> {
  std::vector<DegreeOfFreedom *> Dofs;
  this->getH1Dofs(pointers, Dofs, meshID, order);
  return Dofs;
}

void LinearTriangle::getH1Dofs(PointerCollection &pointers,
                               std::vector<DegreeOfFreedom *> &Dofs,
                               indexType meshID, indexType order) {
  std::vector<DegreeOfFreedom *> tdofs;
  NodeSet *tempSet;
  for (auto i = 0; i < 3; ++i) {
    auto &tempVert = pointers.getGeometryData()->getVertex(this->verts[i]);
    tempSet = tempVert.getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
  if (order > 1) {
    for (auto i = 0; i < 3; ++i) {
      auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[i]);
      tempEdge.getH1DofsInternal(pointers, Dofs, meshID, order);
    }
    this->getH1DofsInternal(pointers, Dofs, meshID, order);
  }
}

void LinearTriangle::getH1DofsInternal(PointerCollection &pointers,
                                       std::vector<DegreeOfFreedom *> &Dofs,
                                       indexType meshID, indexType order) {
  if (order > 2) {
    std::vector<DegreeOfFreedom *> tdofs;
    NodeSet *tempSet;
    tempSet = this->getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
}

auto LinearTriangle::getH1Nodes(PointerCollection &pointers, indexType meshID,
                                indexType order)
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

auto LinearTriangle::getH1NodesInternal(PointerCollection &pointers,
                                        indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  if (order > 2) {
    indexType totnodes = (order - 1) * (order - 2) / 2;
    auto tempnodes = this->getNodesOfSet(pointers, meshID);
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

auto LinearTriangle::getH1Shapes(PointerCollection &pointers, indexType order,
                                 IntegrationPoint &IntegrationPt) -> H1Shapes {

  H1Shapes shapes;
  indexType vertshapes = 3;
  indexType edgeShapes = (order - 1) * 3;
  indexType faceShapes = 0;
  if (order > 2) {
    faceShapes = (order - 1) * (order - 2) / 2;
  }

  indexType numshapes = vertshapes + edgeShapes + faceShapes;
  shapes.shapes.resize(numshapes);
  shapes.shapeDeriv.resize(2, numshapes);

  prec xsi = IntegrationPt.xi;
  prec eta = IntegrationPt.eta;

  prec Lam1 = prec(0.5) * eta + prec(0.5);
  prec Lam2 = -prec(0.5) * xsi - prec(0.5) * eta;
  prec Lam3 = prec(0.5) * xsi + prec(0.5);

  prec dLam1x = prec(0);
  prec dLam2x = -prec(0.5);
  prec dLam3x = prec(0.5);

  prec dLam1e = prec(0.5);
  prec dLam2e = -prec(0.5);
  prec dLam3e = prec(0);

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
      auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[0]);
      prec orientation =
          tempEdge.getEdgeOrientation(this->verts[0], this->verts[1]);

      for (auto i = 0; i < order - 1; ++i) {
        KernelShapes::getShape(tempshape, tempEdgeDerivative,
                               (Lam3 - Lam2) * orientation, i);

        // tempEdgeDerivative *= orientation;

        shapes.shapes(counter) = tempshape * Lam2 * Lam3;
        shapes.shapeDeriv(0, counter) =
            dLam2x * Lam3 * tempshape + Lam2 * dLam3x * tempshape +
            Lam2 * Lam3 * tempEdgeDerivative * (dLam3x - dLam2x) * orientation;
        shapes.shapeDeriv(1, counter) =
            dLam2e * Lam3 * tempshape + Lam2 * dLam3e * tempshape +
            Lam2 * Lam3 * tempEdgeDerivative * (dLam3e - dLam2e) * orientation;
        ++counter;
      }
    }

    // Edge 2
    {
      auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[1]);
      prec orientation =
          tempEdge.getEdgeOrientation(this->verts[1], this->verts[2]);


      for (auto i = 0; i < order - 1; ++i) {
        KernelShapes::getShape(tempshape, tempEdgeDerivative,
                               (Lam1 - Lam3) * orientation, i);

        // tempEdgeDerivative *= orientation;

        shapes.shapes(counter) = tempshape * Lam3 * Lam1;
        shapes.shapeDeriv(0, counter) =
            dLam1x * Lam3 * tempshape + Lam1 * dLam3x * tempshape +
            Lam1 * Lam3 * tempEdgeDerivative * (dLam1x - dLam3x) * orientation;
        shapes.shapeDeriv(1, counter) =
            dLam1e * Lam3 * tempshape + Lam1 * dLam3e * tempshape +
            Lam1 * Lam3 * tempEdgeDerivative * (dLam1e - dLam3e) * orientation;
        ++counter;
      }
    }
    // Edge 3
    {
      auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[2]);
      prec orientation =
          tempEdge.getEdgeOrientation(this->verts[0], this->verts[2]);


      for (auto i = 0; i < order - 1; ++i) {
        KernelShapes::getShape(tempshape, tempEdgeDerivative,
                               (Lam1 - Lam2) * orientation, i);
        
        // tempEdgeDerivative *= orientation;

        shapes.shapes(counter) = tempshape * Lam1 * Lam2;
        shapes.shapeDeriv(0, counter) =
            dLam1x * Lam2 * tempshape + Lam1 * dLam2x * tempshape +
            Lam1 * Lam2 * tempshape * (dLam1x - dLam2x) * orientation;
        shapes.shapeDeriv(1, counter) =
            dLam1e * Lam2 * tempshape + Lam1 * dLam2e * tempshape +
            Lam1 * Lam2 * tempshape * (dLam1e - dLam2e) * orientation;
        ++counter;
      }
    }
  }

  // Face
  if (order > 2) {
    prec tempshape1, tempshape2;
    prec tempShapeDerivative1, tempShapeDerivative2;
    Types::Vector2<prec> Lam_Deriv;

    Lam_Deriv(0) =
        dLam1x * Lam2 * Lam3 + Lam1 * dLam2x * Lam3 + Lam1 * Lam2 * dLam3x;
    Lam_Deriv(1) =
        dLam1e * Lam2 * Lam3 + Lam1 * dLam2e * Lam3 + Lam1 * Lam2 * dLam3e;

    for (auto i = 0; i < order - 2; ++i) {
      for (auto j = 0; j < order - 2 - i; ++j) {
        KernelShapes::getShape(tempshape1, tempShapeDerivative1, Lam3 - Lam2,
                               i);
        KernelShapes::getShape(tempshape2, tempShapeDerivative2, Lam2 - Lam1,
                               j);
        prec lams = Lam1 * Lam2 * Lam3;

        shapes.shapes(counter) = lams * tempshape1 * tempshape2;
        shapes.shapeDeriv(0, counter) =
            Lam_Deriv(0) * tempshape1 * tempshape2 +
            lams * tempShapeDerivative1 * (dLam3x - dLam2x) * tempshape2 +
            lams * tempshape1 * tempShapeDerivative2 * (dLam2x - dLam1x);
        shapes.shapeDeriv(1, counter) =
            Lam_Deriv(1) * tempshape1 * tempshape2 +
            lams * tempShapeDerivative1 * (dLam3e - dLam2e) * tempshape2 +
            lams * tempshape1 * tempShapeDerivative2 * (dLam2e - dLam1e);
        ++counter;
      }
    }
  }

  return shapes;
}

// auto LinearTriangle::getH1ShapesInternal(PointerCollection &pointers,
//                                               indexType order,
//                                               IntegrationPoint
//                                               &IntegrationPt, faceorientation
//                                               orientation)
//     -> H1Shapes {
//   H1Shapes shapes;

//   prec xsi = IntegrationPt.xi;
//   prec eta = IntegrationPt.eta;

//   indexType numShapes = (order - 1) * (order - 2)/2;
//   shapes.shapes.resize(numShapes);
//   shapes.shapeDeriv.resize(2, numShapes);
//   indexType counter = 0;

//   Types::Matrix22<prec> R;

//   if (orientation != faceorientation::p_1) {
//     R.setZero();
//     prec one = prec(1.0);
//     switch (orientation) {
//     case faceorientation::p_2:
//       xsi = -IntegrationPt.eta;
//       eta = IntegrationPt.xi;
//       R(0, 1) = one;
//       R(1, 0) = -one;
//       break;
//     case faceorientation::p_3:
//       xsi = -IntegrationPt.xi;
//       eta = -IntegrationPt.eta;
//       R(0, 0) = -one;
//       R(1, 1) = -one;

//       break;
//     case faceorientation::p_4:
//       xsi = IntegrationPt.eta;
//       eta = -IntegrationPt.xi;
//       R(0, 1) = -one;
//       R(1, 0) = one;
//       break;
//     case faceorientation::n_1:
//       xsi = -IntegrationPt.xi;
//       eta = IntegrationPt.eta;
//       R(0, 0) = -one;
//       R(1, 1) = one;
//       break;
//     case faceorientation::n_2:
//       xsi = -IntegrationPt.eta;
//       eta = -IntegrationPt.xi;
//       R(0, 1) = -one;
//       R(1, 0) = -one;
//       break;
//     case faceorientation::n_3:
//       xsi = IntegrationPt.xi;
//       eta = -IntegrationPt.eta;
//       R(0, 0) = one;
//       R(1, 1) = -one;
//       break;
//     case faceorientation::n_4:
//       xsi = IntegrationPt.eta;
//       eta = IntegrationPt.xi;
//       R(0, 1) = one;
//       R(1, 0) = one;
//       break;
//     default: {
//     }
//     }
//   }

//   for (auto i = 2; i <= order; i++) {
//     auto xiShape = LobattoShapes::getShape(xsi, i);
//     // xiShape.shapeDerivative *= factXi;
//     for (auto j = 2; j <= order; j++) {
//       auto etaShape = LobattoShapes::getShape(eta, j);
//       // etaShape.shapeDerivative *= factEta;
//       shapes.shapes(counter) = xiShape.shapeValue * etaShape.shapeValue;
//       shapes.shapeDeriv(0, counter) =
//           xiShape.shapeDerivative * etaShape.shapeValue;
//       shapes.shapeDeriv(1, counter) =
//           xiShape.shapeValue * etaShape.shapeDerivative;
//       ++counter;
//     }
//   }
//   if (orientation != faceorientation::p_1) {
//     shapes.shapeDeriv = R * shapes.shapeDeriv;
//   }

//   return shapes;
// }

// ----------------HDiv------------------

auto LinearTriangle::getHDivNodes(PointerCollection &pointers, indexType meshID,
                                  indexType order)
    -> std::vector<GenericNodes *> {
  return {};
}

auto LinearTriangle::getHDivNodesInternal(PointerCollection &pointers,
                                          indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  if (order > 1) {
    indexType totnodes = (order - 1) * (order - 1);
    auto tempnodes = this->getNodesOfSet(pointers, meshID);
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

void LinearTriangle::setHDivShapes(PointerCollection &pointers,
                                   indexType meshId, indexType order,
                                   NodeTypes type) {
  indexType numNodes = order + 1;
  for (auto i = 0; i < 3; ++i) {
    auto &edgeTemp = pointers.getGeometryData()->getEdge(this->edges[i]);
    edgeTemp.setNodeSet(pointers, meshId, numNodes, type);
  }
  if (order > 1) {
    numNodes = order * order - 1;
    this->setNodeSet(pointers, meshId, numNodes, type);
  }
}

void LinearTriangle::getHDivDofs(PointerCollection &pointers,
                                 std::vector<DegreeOfFreedom *> &Dofs,
                                 indexType meshID, indexType order,
                                 NodeTypes type) {
  std::vector<DegreeOfFreedom *> tdofs;
  NodeSet *tempSet;
  for (auto i = 0; i < 3; ++i) {
    auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[i]);
    tempSet = tempEdge.getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
  if (order > 1) {
    tempSet = this->getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
}

auto LinearTriangle::getHDivShapes(PointerCollection &pointers, indexType order,
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
    auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[0]);
    prec orientation =
        tempEdge.getEdgeOrientation(this->verts[0], this->verts[1]);
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
    auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[1]);
    prec orientation =
        tempEdge.getEdgeOrientation(this->verts[1], this->verts[2]);
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
    auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[2]);
    prec orientation =
        tempEdge.getEdgeOrientation(this->verts[2], this->verts[0]);
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

void LinearTriangle::setAllNodeBoundaryConditionMeshId(
    PointerCollection &pointers, indexType meshId, indexType dof) {
  Base::setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);

  for (auto i = 0; i < 3; i++) {
    auto &tempEdge = pointers.getGeometryData()->getEdge(this->edges[i]);
    tempEdge.setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);
  }
}

void LinearTriangle::setLoad(PointerCollection &pointers, indexType meshid,
                             ShapeFunctionTypes shapeType, indexType shapeOrder,
                             Types::VectorX<prec> &Loads, indexType propNumber,
                             Types::VectorX<prec> &direction, bool local,
                             bool add) {

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
      for (auto j = 0; j < 3; ++j) {
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

void LinearTriangle::setBoundaryCondition(PointerCollection &pointers,
                                          indexType meshId, indexType order,
                                          ShapeFunctionTypes shapeType,
                                          Types::Vector3<indexType> &dofs,
                                          bool set) {
  if (shapeType == ShapeFunctionTypes::H1) {
    auto nodes = this->getH1Nodes(pointers, meshId, order);
    for (auto i : nodes) {
      if (set) {
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            i->setBoundaryCondition(pointers, j);
          } else {
            i->unsetBoundaryCondition(pointers, dofs(j));
          }
        }
      } else {
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            i->setBoundaryCondition(pointers, dofs(j));
          }
        }
      }
    }
  }
}

void LinearTriangle::geometryToParaview(PointerCollection &pointers,
                                        vtkPlotInterface &paraviewAdapter,
                                        indexType mainMesh, indexType subMesh) {
  indexType numPoints = 3;
  std::vector<indexType> points(numPoints);
  points.clear();
  for (auto i : this->verts) {
    auto &vert = pointers.getGeometryData()->getVertex(i);
    vert.geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
    points.push_back(i);
  }

  paraviewAdapter.addCell(mainMesh, subMesh, this->id, 1, points, 3,
                          VTK_TRIANGLE);
}

void LinearTriangle::computeWeightsParaview(PointerCollection &pointers,
                                            vtkPlotInterface &paraviewAdapter,
                                            indexType mainMesh,
                                            indexType subMesh) {

  auto GP = this->getIntegrationPoints(pointers,-1);
  GP.setOrder(2);

  for (auto i : GP) {
    auto jaco = this->getJacobian(pointers, i);
    auto shapes = this->getH1Shapes(pointers, 1, i);
    prec dA = jaco.determinant() * i.weight;

    for (auto i = 0; i < 3; ++i) {
      auto &vert = pointers.getGeometryData()->getVertex(this->verts[i]);
      std::vector<prec> val;
      val.push_back(shapes.shapes(i) * dA);
      paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, val, vert.getId(),
                                           1, paraviewNames::weightName());
    }
  }
}

void LinearTriangle::H1SolutionToParaview(PointerCollection &pointers,
                                          vtkPlotInterface &paraviewAdapter,
                                          indexType mainMesh, indexType subMesh,
                                          indexType meshId, indexType order,
                                          std::string &name) {
  std::vector<DegreeOfFreedom *> Dofs;
  this->getH1Dofs(pointers, Dofs, meshId, order);
  auto solution = pointers.getSolutionState()->getSolution(Dofs);
  for (auto i = 0; i < 3; ++i) {
    auto &V = pointers.getGeometryData()->getVertex(this->verts[i]);
    std::vector<prec> sol(3);
    for (auto j = 0; j < 3; ++j) {
      sol[j] = solution(3 * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, V.getId(), sol, 3, name);
  }
}

void LinearTriangle::H1DataToParaview(PointerCollection &pointers,
                                      vtkPlotInterface &paraviewAdapter,
                                      indexType mainMesh, indexType subMesh,
                                      Types::VectorX<prec> &Data,
                                      indexType numberComponents,
                                      indexType order, std::string &name) {
  for (auto i = 0; i < 3; ++i) {
    auto &V = pointers.getGeometryData()->getVertex(this->verts[i]);
    std::vector<prec> sol(numberComponents);
    for (auto j = 0; j < numberComponents; ++j) {
      sol[j] = Data(numberComponents * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, V.getId(), sol,
                                 numberComponents, name);
  }
}

void LinearTriangle::projectDataToParaviewVertices(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {

  auto shapes = this->getH1Shapes(pointers, 1, IntegrationPt);
  std::vector<prec> vals(numberComponents);

  auto jaco = this->getJacobian(pointers, IntegrationPt);
  auto dA = jaco.determinant() * IntegrationPt.weight;
  for (auto i = 0; i < 3; ++i) {
    auto &V = pointers.getGeometryData()->getVertex(this->verts[i]);
    for (auto j = 0; j < numberComponents; ++j) {
      vals[j] = data(j) * shapes.shapes(i) * dA;
    }
    paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, vals, V.getId(),
                                         numberComponents, name);
  }
}

void LinearTriangle::flip() {
  throw std::runtime_error(
      "ERROR in LinearTriangle::flip: Method not implemented!");
}

void LinearTriangle::rotate(indexType n) {
  throw std::runtime_error(
      "ERROR in LinearTriangle::rotate: Method not implemented!");
}

auto LinearTriangle::computeMeanCoordinate(PointerCollection &pointers)
    -> Types::Vector3<prec> {
  throw std::runtime_error(
      "ERROR in LinearTriangle::computeMeanCoordinate: "
      "Method not implemented!");
  return Types::Vector3<prec>();
}

const GeometryTypes LinearTriangle::type = GeometryTypes::LinearTriangle;

} // namespace HierAMuS::Geometry
