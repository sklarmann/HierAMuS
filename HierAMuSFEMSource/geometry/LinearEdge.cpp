// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "geometry/Base.h"
#include "geometry/GeometryTypes.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

#include "geometry/GeometryData.h"
#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <equations/NodeSet.h>
#include <geometry/GeometryData.h>
#include <geometry/LinearEdge.h>
#include <geometry/Vertex.h>
#include <pointercollection/pointercollection.h>

#include <equations/DegreeOfFreedom.h>

#include <loads/LoadList.h>

#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

#include <shapefunctions/LagrangeShape.h>
#include <stdexcept>
#include <types/MatrixTypes.h>
#include <vtkCellType.h>

#include <iomanip>
#include <sstream>

namespace HierAMuS::Geometry {
using std::vector;

LinearEdge::LinearEdge() : m_verts({-1, -1}) {}

LinearEdge::~LinearEdge() = default;

auto LinearEdge::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::LinearEdge::type;
}

void LinearEdge::getVerts(std::vector<indexType> &vertsOut) {
  vertsOut.resize(2);
  for (auto i = 0; i < 2; ++i) {
    vertsOut[i] = this->m_verts[i];
  }
}

void LinearEdge::getVerts(PointerCollection &pointers,
                          std::vector<Base *> &vertsOut) {
  vertsOut.clear();
  for (auto i = 0; i < 2; ++i) {
    vertsOut.push_back(
        &pointers.getGeometryData()->getVertex(this->m_verts[i]));
  }
}

auto LinearEdge::getVertex(PointerCollection &pointers, indexType number)
    -> Vertex & {
  if (number > 2) {
    throw std::runtime_error(
        "Requested number of vertex too large. Element has only two vertices!");
  }
  return pointers.getGeometryData()->getVertex(this->m_verts[number]);
}

auto LinearEdge::getVerts() -> std::vector<indexType> {
  return std::vector<indexType>(m_verts.begin(), m_verts.end());
}

auto LinearEdge::getVertexNumber(indexType number) -> indexType {
  return this->m_verts[number];
}

auto LinearEdge::hasVertices(indexType v1, indexType v2) -> bool {
  return (this->m_verts[0] == v1 && this->m_verts[1] == v2) ||
         (this->m_verts[0] == v2 && this->m_verts[1] == v1);
}

void LinearEdge::setVerts(GeometryData &geoData,
                          std::vector<indexType> &vertsIn) {
  if (vertsIn.size() == 2) {
    for (auto i = 0; i < 2; ++i) {
      this->m_verts[i] = vertsIn[i];
      auto &Vert = geoData.getVertex(vertsIn[i]);
      Vert.connectEdge(this->getId());
    }
  }
}

void LinearEdge::getNodes(PointerCollection &pointers,
                          std::vector<GenericNodes *> &nodeVector,
                          indexType meshId) {
  nodeVector.clear();
  std::vector<GenericNodes *> tempNodes;
  for (auto &i : this->m_verts) {
    auto &tempElem = pointers.getGeometryData()->getVertex(i);
    tempElem.getNodesOfSet(pointers, tempNodes, meshId);
    nodeVector.insert(nodeVector.end(), tempNodes.begin(), tempNodes.end());
  }
  this->getNodesOfSet(pointers, tempNodes, meshId);
  nodeVector.insert(nodeVector.end(), tempNodes.begin(), tempNodes.end());
}

inline void LinearEdge::print(PointerCollection &pointers) {
  auto &Logger = pointers.getSPDLogger();

  Logger.debug("Linear Edge id: {:d}", this->id);
  Logger.debug("Vertices: {:d}, {:d}", m_verts[0], m_verts[1]);

  this->printEqInfo(pointers);
}

auto LinearEdge::getEdgeOrientation(indexType startVertex, indexType endVertex)
    -> prec {
  prec val;
  startVertex == this->m_verts[0] ? val = prec(1.0) : val = prec(-1.0);
  // val = prec(1.0);
  return val;
}

auto LinearEdge::getCoordinates(PointerCollection &pointers, prec xi)
    -> Types::Vector3<prec> {
  prec v1 = prec(1) - xi;
  v1 /= prec(2);
  prec v2 = prec(1) + xi;
  v2 /= prec(2);

  Types::Vector3<prec> c1;
  Types::Vector3<prec> c2;
  Types::Vector3<prec> c3;

  c1 = pointers.getGeometryData()->getVertex(this->m_verts[0]).getCoordinates();
  c2 = pointers.getGeometryData()->getVertex(this->m_verts[1]).getCoordinates();

  c3 = v1 * c1 + v2 * c2;

  return c3;
}

auto LinearEdge::getCoordinates(PointerCollection &pointers,
                                IntegrationPoint &IntPoint)
    -> Types::Vector3<prec> {
  prec v1 = prec(1) - IntPoint.xi;
  v1 /= prec(2);
  prec v2 = prec(1) + IntPoint.xi;
  v2 /= prec(2);

  Types::Vector3<prec> c1;
  Types::Vector3<prec> c2;
  Types::Vector3<prec> c3;

  c1 = pointers.getGeometryData()->getVertex(this->m_verts[0]).getCoordinates();
  c2 = pointers.getGeometryData()->getVertex(this->m_verts[1]).getCoordinates();

  c3 = v1 * c1 + v2 * c2;

  return c3;
}

auto LinearEdge::getIntegrationPoints(PointerCollection &pointers,
                                      indexType elementId)
    -> IntegrationPoints {

  auto intpoints = HierAMuS::PointerCollection::getIntegrationPoints(elementId);
  intpoints.setType(IntegrationType::Gauss1D);
  return intpoints;
}

void LinearEdge::setH1Shapes(PointerCollection &pointers, indexType meshId,
                             indexType order, NodeTypes type) {

  for (auto i = 0; i < 2; ++i) {
    auto &tempGeo = pointers.getGeometryData()->getVertex(this->m_verts[i]);
    tempGeo.setNodeSet(pointers, meshId, 1, type);
  }
  this->setH1ShapesInternal(pointers, meshId, order, type);
}

void LinearEdge::setH1ShapesInternal(PointerCollection &pointers,
                                     indexType meshId, indexType order,
                                     NodeTypes type) {
  if (order > 1) {
    this->setNodeSet(pointers, meshId, order - 1, type);
  }
}

void LinearEdge::getH1Dofs(PointerCollection &pointers,
                           std::vector<DegreeOfFreedom *> &Dofs,
                           indexType meshID, indexType order) {
  std::vector<DegreeOfFreedom *> tdofs;
  NodeSet *tempSet;
  for (auto i = 0; i < 2; ++i) {
    auto &tempVert = pointers.getGeometryData()->getVertex(this->m_verts[i]);
    tempSet = tempVert.getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
  this->getH1DofsInternal(pointers, Dofs, meshID, order);
}

void LinearEdge::getH1DofsInternal(PointerCollection &pointers,
                                   std::vector<DegreeOfFreedom *> &Dofs,
                                   indexType meshID, indexType order) {
  if (order > 1) {
    std::vector<DegreeOfFreedom *> tdofs;
    NodeSet *tempSet;
    tempSet = this->getSetMeshId(pointers, meshID);
    tdofs = tempSet->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), tdofs.begin(), tdofs.end());
  }
}

auto LinearEdge::getH1Nodes(PointerCollection &pointers, indexType meshID,
                            indexType order) -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  for (auto i : this->m_verts) {
    auto &tempVert = pointers.getGeometryData()->getVertex(i);
    auto tempNodes = tempVert.getNodesOfSet(pointers, meshID);
    nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
  }
  if (order > 1) {
    auto tempnodes = this->getH1NodesInternal(pointers, meshID, order);
    nodes.insert(nodes.end(), tempnodes.begin(), tempnodes.end());
  }
  return nodes;
}

auto LinearEdge::getH1NodesInternal(PointerCollection &pointers,
                                    indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  if (order > 1) {
    auto tempnodes = this->getNodesOfSet(pointers, meshID);
    nodes.insert(nodes.end(), tempnodes.begin(), tempnodes.end());
    if (tempnodes.size() != order - 1) {
      std::stringstream ss;
      ss << "Error in Linear Edge geometry element with id: " << this->id
         << " and order: " << order << " has " << tempnodes.size()
         << " nodes instead of " << order + 1 << " nodes.";
      throw std::runtime_error(ss.str());
    }
  }
  return nodes;
}

/**
 * @brief Get local H1 shape functions of Edge.
 *
 * @param pointers Pointers to global data.
 * @param order Order of the shape functions. 1: Linear, 2: Quadratic ....
 * @param shape Shape functions.
 * @param shapeDerivative Local derivatives of the shape functions.
 * @param xsi Local parameter xsi to evaluate the shape functions and
 * derivatives.
 */

void LinearEdge::getH1Shapes(PointerCollection &pointers, indexType order,
                             Types::VectorX<prec> &shape,
                             Types::VectorX<prec> &shapeDerivative, prec xsi) {
  shape.resize(order + 1);
  shapeDerivative.resize(order + 1);
  for (auto i = 0; i < 2; ++i) {
    // LobattoShape(shape(i), shapeDerivative(i), xsi, i);
    LobattoShapes::getShape(shape(i), shapeDerivative(i), xsi, i);
  }
  Types::VectorX<prec> tempshape;
  Types::VectorX<prec> tempshapeDerivative;
  this->getH1ShapesInternal(pointers, order, tempshape, tempshapeDerivative,
                            xsi);
  for (auto i = 0; i < tempshape.rows(); ++i) {
    shape(i + 2) = tempshape(i);
    shapeDerivative(i + 2) = tempshapeDerivative(i);
  }
}

void LinearEdge::getH1ShapesInternal(PointerCollection &pointers,
                                     indexType order,
                                     Types::VectorX<prec> &shape,
                                     Types::VectorX<prec> &shapeDerivative,
                                     prec xsi) {
  shape.resize(order - 1);
  shapeDerivative.resize(order - 1);
  for (auto i = 0; i < order - 1; ++i) {
    // LobattoShape(shape(i), shapeDerivative(i), xsi, i + 2);
    LobattoShapes::getShape(shape(i), shapeDerivative(i), xsi, i + 2);
  }
}

auto LinearEdge::getH1Shapes(PointerCollection &pointers, indexType order,
                             IntegrationPoint &integration_point) -> H1Shapes {
  H1Shapes shapes;
  shapes.shapes.resize(order + 1);
  shapes.shapeDeriv.resize(1, order + 1);

  for (auto i = 0; i < 2; ++i) {
    // LobattoShape(shape(i), shapeDerivative(i), xsi, i);
    LobattoShapes::getShape(shapes.shapes(i), shapes.shapeDeriv(0, i),
                            integration_point.xi, i);
  }
  if (order > 1) {
    auto addShapes =
        this->getH1ShapesInternal(pointers, order, integration_point);
    shapes.shapes.tail(order - 1) = addShapes.shapes;
    // shapes.shapes.block(0, 2, 1, order - 1) =
    //     addShapes.shapes.block(0, 0, 1, order - 1);
    shapes.shapeDeriv.block(0, 2, 1, order - 1) =
        addShapes.shapeDeriv.block(0, 0, 1, order - 1);
  }

  return shapes;
}

auto LinearEdge::getH1ShapesInternal(PointerCollection &pointers,
                                     indexType order,
                                     IntegrationPoint &integration_point)
    -> H1Shapes {

  H1Shapes shapes;
  if (order > 1) {
    shapes.shapes.resize(order - 1);
    shapes.shapeDeriv.resize(1, order - 1);
    for (auto i = 2; i < order + 1; ++i) {
      HierAMuS::LobattoShapes::getShape(shapes.shapes(i - 2),
                                        shapes.shapeDeriv(0, i - 2),
                                        integration_point.xi, i);
    }
  }

  return shapes;
}

auto LinearEdge::getJacobian(PointerCollection &pointers, prec xi) -> prec {

  prec jac;

  Types::Vector3<prec> dirVec;
  dirVec =
      pointers.getGeometryData()->getVertex(this->m_verts[1]).getCoordinates();
  dirVec -=
      pointers.getGeometryData()->getVertex(this->m_verts[0]).getCoordinates();

  jac = dirVec.norm() / prec(2);

  return jac;
}

auto LinearEdge::getJacobian(PointerCollection &pointers,
                             IntegrationPoint &IntegrationPt)
    -> Types::MatrixXX<prec> {

  Types::MatrixXX<prec> jacob;
  jacob.resize(1, 1);
  Types::Vector3<prec> dirVec;
  dirVec = pointers.getGeometryData()->getVertex(m_verts[1]).getCoordinates();
  dirVec -= pointers.getGeometryData()->getVertex(m_verts[0]).getCoordinates();

  jacob(0, 0) = dirVec.norm() / static_cast<prec>(2);

  return jacob;
}

void LinearEdge::setLoad(PointerCollection &pointers, indexType meshid,
                         ShapeFunctionTypes shapeType, indexType shapeOrder,
                         Types::VectorX<prec> &Loads, indexType propNumber,
                         Types::VectorX<prec> &direction, bool local,
                         bool add) {
  if (shapeType == ShapeFunctionTypes::HDiv ||
      shapeType == ShapeFunctionTypes::HCurl) {
    std::vector<DegreeOfFreedom *> tempDofs;
    GenericNodes *tnode;
    std::vector<GenericNodes *> tempNodes;

    this->getNodes(pointers, tempNodes, meshid);

    tnode = tempNodes[0];
    tnode->getDegreesOfFreedom(pointers, tempDofs);
    for (auto i = 0; i < 3; ++i) {
      pointers.getLoadList()->setLoad(propNumber, tempDofs[i]->getId(),
                                      Loads(i), false);
    }

  } else {
    auto GP = this->getIntegrationPoints(pointers, -1);
    GP.setOrder(shapeOrder);

    Types::VectorX<prec> LoadShape;
    Types::VectorX<prec> LoadShapederiv;

    indexType loadOrder;
    loadOrder = Loads.size() / 3 - 1;

    if (loadOrder == 0) {
      LoadShape.resize(1);
      LoadShape(0) = prec(1);
    }

    std::vector<GenericNodes *> tempNodes;
    this->getNodes(pointers, tempNodes, meshid);

    if (shapeOrder > tempNodes.size() - 1) {
      shapeOrder = tempNodes.size() - 1;
    }

    std::vector<DegreeOfFreedom *> tempDofs;
    GenericNodes *tnode;

    if (!add) {
      for (auto &tempNode : tempNodes) {
        tnode = tempNode;
        tnode->getDegreesOfFreedom(pointers, tempDofs);
        for (auto i = 0; i < 3; ++i) {
          pointers.getLoadList()->setLoad(propNumber, tempDofs[i]->getId(),
                                          prec(0), false);
        }
      }
    }

    for (auto i : GP) {
      auto jacobi = this->getJacobian(pointers, i);
      auto h1Shapes = this->getH1Shapes(pointers, shapeOrder, i);
      prec dA = jacobi.determinant() * i.weight;

      Types::Vector3<prec> localLoad;
      if (Loads.size() / 3 > 1) {
        localLoad.setZero();
        indexType nLoadShapes;
        nLoadShapes = Loads.size() / 3;
        for (indexType nl = 0; nl < nLoadShapes; ++nl) {
          prec sh = prec(0);
          prec dsh = prec(0);
          LagrangeShape(sh, dsh, i.xi, nLoadShapes - 1, nl);
          localLoad += sh * localLoad.block(3 * nl, 0, 3, 1);
        }
      } else {
        localLoad = Loads;
      }
      for (auto nn = 0; nn < shapeOrder + 1; ++nn) {
        tnode = tempNodes[nn];

        tnode->getDegreesOfFreedom(pointers, tempDofs);
        for (auto j = 0; j < 3; ++j) {
          prec loadValue = localLoad(j);
          loadValue *= h1Shapes.shapes(nn);
          loadValue *= dA;
          pointers.getLoadList()->setLoad(propNumber, tempDofs[j]->getId(),
                                          loadValue, true);
        }
      }
    }
  }
}

void LinearEdge::setPrescribedSolution(
    PointerCollection &pointers, indexType meshid, ShapeFunctionTypes shapeType,
    indexType shapeOrder, Types::VectorX<prec> &Solution, indexType propNumber,
    Types::VectorX<prec> &direction, bool local, bool add) {
  if (shapeType == ShapeFunctionTypes::H1) {
    std::vector<GenericNodes *> tempNodes =
        this->getH1Nodes(pointers, meshid, shapeOrder);
    for (indexType i = 0; i < tempNodes.size(); i++) {
      for (indexType j = 0; j < 3; ++j) {
        if (abs(Solution(j)) > std::numeric_limits<prec>::epsilon() * 1000) {
          auto &dof = tempNodes[i]->getDegreeOfFreedom(j);
          dof.setStatus(dofStatus::inactive);
          if (j < 2)
            pointers.getPrescribedDisplacements()->setLoad(
                propNumber, dof.getId(), Solution(j), add);
        }
      }
    }
  } else {
    pointers.getSPDLogger().warn(
        "In LinearEdge::setPrescribedSolution: Trying to set prescibed "
        "solution on Edge {} with an unsupported shape function type.",
        this->id);
  }
}

void LinearEdge::geometryToParaview(PointerCollection &pointers,
                                    vtkPlotInterface &paraviewAdapter,
                                    indexType mainMesh, indexType subMesh) {
  std::vector<indexType> points(2);
  points = {m_verts[0], m_verts[1]};
  for (auto i : this->m_verts) {
    auto &vert = pointers.getGeometryData()->getVertex(i);
    vert.geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
  }
  paraviewAdapter.addCell(mainMesh, subMesh, this->id, 1, points, 2, VTK_LINE);
}

void LinearEdge::computeWeightsParaview(PointerCollection &pointers,
                                        vtkPlotInterface &paraviewAdapter,
                                        indexType mainMesh, indexType subMesh) {
  auto GP = this->getIntegrationPoints(pointers, -1);
  GP.setOrder(2);
  for (auto i : GP) {
    auto jaco = this->getJacobian(pointers, i);
    auto shapes = this->getH1Shapes(pointers, 1, i);
    prec dA = jaco.determinant() * i.weight;
    for (auto j = 0; j < 2; ++j) {
      auto &vert = pointers.getGeometryData()->getVertex(this->m_verts[j]);
      std::vector<prec> val;
      val.push_back(shapes.shapes(j) * dA);
      paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, val, vert.getId(),
                                           1, paraviewNames::weightName());
    }
  }
}

void LinearEdge::H1SolutionToParaview(PointerCollection &pointers,
                                      vtkPlotInterface &paraviewAdapter,
                                      indexType mainMesh, indexType subMesh,
                                      indexType meshId, indexType order,
                                      std::string &name) {

  std::vector<DegreeOfFreedom *> Dofs;
  this->getH1Dofs(pointers, Dofs, meshId, order);
  auto solution = pointers.getSolutionState()->getSolution(Dofs);
  for (auto i = 0; i < 2; ++i) {
    auto &V = pointers.getGeometryData()->getVertex(this->m_verts[i]);
    std::vector<prec> sol(3);
    for (auto j = 0; j < 3; ++j) {
      sol[j] = solution(3 * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, V.getId(), sol, 3, name);
  }
}

void LinearEdge::H1DataToParaview(PointerCollection &pointers,
                                  vtkPlotInterface &paraviewAdapter,
                                  indexType mainMesh, indexType subMesh,
                                  Types::VectorX<prec> &Data,
                                  indexType numberComponents, indexType order,
                                  std::string &name) {
  for (auto i = 0; i < 2; ++i) {
    auto &V = pointers.getGeometryData()->getVertex(this->m_verts[i]);
    std::vector<prec> sol(numberComponents);
    for (auto j = 0; j < numberComponents; ++j) {
      sol[j] = Data(numberComponents * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, V.getId(), sol,
                                 numberComponents, name);
  }
}

void LinearEdge::projectDataToParaviewVertices(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {
  auto shapes = this->getH1Shapes(pointers, 1, IntegrationPt);
  std::vector<prec> vals(numberComponents);

  auto jaco = this->getJacobian(pointers, IntegrationPt);
  auto dA = jaco.determinant() * IntegrationPt.weight;
  for (auto i = 0; i < 2; ++i) {
    auto &V = pointers.getGeometryData()->getVertex(this->m_verts[i]);
    for (auto j = 0; j < numberComponents; ++j) {
      vals[j] = data(j) * shapes.shapes(i) * dA;
    }
    paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, vals, V.getId(),
                                         numberComponents, name);
  }
}

auto LinearEdge::getDirectionVector(PointerCollection &pointers)
    -> Types::Vector3<prec> {
  auto &V1 = pointers.getGeometryData()->getVertex(this->m_verts[0]);
  auto &V2 = pointers.getGeometryData()->getVertex(this->m_verts[1]);
  Types::Vector3<prec> dir;
  dir = V2.getCoordinates() - V1.getCoordinates();
  return dir;
}

auto LinearEdge::getA1Vector(PointerCollection &pointers,
                             IntegrationPoint &integration_point)
    -> Types::Vector3<prec> {

  Types::Vector3<prec> dir;
  auto &V1 = pointers.getGeometryData()->getVertex(m_verts[0]);
  auto &V2 = pointers.getGeometryData()->getVertex(m_verts[1]);
  dir = V2.getCoordinates() - V1.getCoordinates();
  dir /= dir.norm();
  return dir;
}

void LinearEdge::setBoundaryCondition(PointerCollection &pointers,
                                      indexType meshId, indexType order,
                                      ShapeFunctionTypes shapeType,
                                      Types::Vector3<indexType> &dofs,
                                      bool set) {

  if (shapeType == ShapeFunctionTypes::H1) {
    std::vector<GenericNodes *> tnodes;
    this->getNodes(pointers, tnodes, meshId);
    indexType nnodes = order + 1;
    tnodes.size() < nnodes ? nnodes = tnodes.size() : nnodes = nnodes;
    if (set) {
      for (auto i = 0; i < nnodes; ++i) {
        GenericNodes *tnode;
        tnode = tnodes[i];
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            tnode->setBoundaryCondition(pointers, j);
          } else {
            tnode->unsetBoundaryCondition(pointers, j);
          }
        }
      }
    } else {
      for (auto i = 0; i < nnodes; ++i) {
        GenericNodes *tnode;
        tnode = tnodes[i];
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            tnode->setBoundaryCondition(pointers, j);
          }
        }
      }
    }

  } else if (shapeType == ShapeFunctionTypes::H0 ||
             shapeType == ShapeFunctionTypes::HCurl ||
             shapeType == ShapeFunctionTypes::HCurl) {
    std::vector<GenericNodes *> tnodes;
    this->getNodesOfSet(pointers, tnodes, meshId);
    indexType nnodes = order - 1;
    tnodes.size() < nnodes ? nnodes = tnodes.size() : nnodes = nnodes;
    if (set) {
      for (auto i = 0; i < nnodes; ++i) {
        GenericNodes *tnode;
        tnode = tnodes[i];
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            tnode->setBoundaryCondition(pointers, j);
          } else {
            tnode->unsetBoundaryCondition(pointers, j);
          }
        }
      }
    } else {
      for (auto i = 0; i < nnodes; ++i) {
        GenericNodes *tnode;
        tnode = tnodes[i];
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            tnode->setBoundaryCondition(pointers, j);
          }
        }
      }
    }
  }
}

void LinearEdge::getNodesInternal(PointerCollection &pointers,
                                  std::vector<GenericNodes *> &nodeVector,
                                  indexType meshId) {
  nodeVector.clear();
  this->getNodesOfSet(pointers, nodeVector, meshId);
}

const GeometryTypes LinearEdge::type = GeometryTypes::LinearEdge;

void LinearEdge::setAllNodeBoundaryConditionMeshId(PointerCollection &pointers,
                                                   indexType meshId,
                                                   indexType dof) {
  Base::setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);

  for (auto i : this->m_verts) {
    auto &V = pointers.getGeometryData()->getVertex(i);
    V.setAllNodeBoundaryConditionMeshId(pointers, meshId, dof);
  }
}

void LinearEdge::flip() {
  indexType temp = this->m_verts[0];
  this->m_verts[0] = this->m_verts[1];
  this->m_verts[1] = temp;
}

} // namespace HierAMuS::Geometry
