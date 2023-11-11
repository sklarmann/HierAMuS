// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "geometry/GeometryBaseData.h"
#include "geometry/GeometryTypes.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "shapefunctions/LobattoShapes.h"

#include "geometry/GeometryData.h"
#include "geometry/Edges/LinearEdgeData.h"
#include "geometry/Edges/LinearEdgeRuntime.h"
#include "geometry/VertexData.h"

#include "plot/vtkplotClass.h"


#include "LoadList.h"

#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

#include <shapefunctions/LagrangeShape.h>
#include <stdexcept>
#include <types/MatrixTypes.h>
#include <vtkCellType.h>

#include <iomanip>
#include <sstream>

#include "shapefunctions/IntegrationsPoints/helperClasses/IntegrationPointsManagement.h"

//Equations
#include "GenericNodes.h"


namespace HierAMuS::Geometry {
using std::vector;

LinearEdgeData::LinearEdgeData() : EdgesDataInterface() {}

LinearEdgeData::~LinearEdgeData() = default;

auto LinearEdgeData::getRuntimeObject(GeometryData &geoData)
    -> std::shared_ptr<EdgesRuntime> {
  return std::make_shared<LinearEdgeRuntime>(geoData, *this);
}

auto LinearEdgeData::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::LinearEdgeData::type;
}



void LinearEdgeData::getNodes(std::vector<GenericNodes *> &nodeVector,
                              indexType meshId) {
  nodeVector.clear();
  for (auto &i : this->m_verts_pointers) {
    auto tempNodes = i->getNodesOfSet(meshId);
    nodeVector.insert(nodeVector.end(), tempNodes.begin(), tempNodes.end());
  }
  auto tempNodes = this->getNodesOfSet(meshId);
  nodeVector.insert(nodeVector.end(), tempNodes.begin(), tempNodes.end());
}

auto LinearEdgeData::getH1NodesList(indexType meshID, indexType order)
    -> MeshIdNodeList {
  MeshIdNodeList tempList(meshID);
  tempList.reserve(3);
  for (auto &i : this->m_verts_pointers) {
    auto ll = i->getNodeSetNodeListMeshId(meshID);
    tempList.add(ll);
  }
  if (order > 1) {
    auto ll = this->getNodeSetNodeListMeshId(meshID);
    tempList.add(ll);
  }
  return tempList;
}


auto LinearEdgeData::getEdgeOrientation(indexType startVertex,
                                        indexType endVertex) -> prec {
  prec val;
  startVertex == this->m_verts[0] ? val = prec(1.0) : val = prec(-1.0);
  // val = prec(1.0);
  return val;
}

auto LinearEdgeData::getCoordinates(prec xi)
    -> Types::Vector3<prec> {
  prec v1 = prec(1) - xi;
  v1 /= prec(2);
  prec v2 = prec(1) + xi;
  v2 /= prec(2);

  Types::Vector3<prec> c1;
  Types::Vector3<prec> c2;
  Types::Vector3<prec> c3;

  c1 = m_verts_pointers[0]->getCoordinates();
  c2 = m_verts_pointers[1]->getCoordinates();

  c3 = v1 * c1 + v2 * c2;

  return c3;
}

auto LinearEdgeData::getCoordinates(IntegrationPoint &IntPoint)
    -> Types::Vector3<prec> {
  prec v1 = prec(1) - IntPoint.xi;
  v1 /= prec(2);
  prec v2 = prec(1) + IntPoint.xi;
  v2 /= prec(2);

  Types::Vector3<prec> c1;
  Types::Vector3<prec> c2;
  Types::Vector3<prec> c3;

  c1 = m_verts_pointers[0]->getCoordinates();
  c2 = m_verts_pointers[1]->getCoordinates();

  c3 = v1 * c1 + v2 * c2;

  return c3;
}

auto LinearEdgeData::getIntegrationPoints(indexType elementId) -> IntegrationPoints {

  auto intpoints = IntegrationPointsManagement::getIntegrationsPoints(elementId);
  intpoints.setType(IntegrationType::Gauss1D);
  return intpoints;
}

void LinearEdgeData::setH1Shapes(indexType meshId,
                                 indexType order, NodeTypes type) {

  for (auto i = 0; i < 2; ++i) {
    m_verts_pointers[i]->setNodeSet(meshId, 1, type);
  }
  this->setH1ShapesInternal(meshId, order, type);
}

void LinearEdgeData::setH1ShapesInternal(indexType meshId, indexType order,
                                         NodeTypes type) {
  if (order > 1) {
    this->setNodeSet(meshId, order - 1, type);
  }
}

void LinearEdgeData::getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs,
                               indexType meshID, indexType order) {

  auto ll = this->getH1NodesList(meshID, order);
  Dofs = ll.getDegreesOfFreedom();
}

void LinearEdgeData::getH1DofsInternal(std::vector<DegreeOfFreedom *> &Dofs,
                                       indexType meshID, indexType order) {
  if (order > 1) {
    auto set = this->getNodeSetNodeListMeshId(meshID);
    set.addDofsToVector(Dofs);
  }
}

auto LinearEdgeData::getH1Nodes(indexType meshID,
                                indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  for (auto &i : this->m_verts_pointers) {
    auto tempNodes = i->getNodesOfSet(meshID);
    nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
  }
  if (order > 1) {
    auto tempnodes = this->getH1NodesInternal(meshID, order);
    nodes.insert(nodes.end(), tempnodes.begin(), tempnodes.end());
  }
  return nodes;
}

auto LinearEdgeData::getH1NodesInternal(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  if (order > 1) {
    auto tempnodes = this->getNodesOfSet(meshID);
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

void LinearEdgeData::getH1Shapes(indexType order,
                                 Types::VectorX<prec> &shape,
                                 Types::VectorX<prec> &shapeDerivative,
                                 prec xsi) {
  shape.resize(order + 1);
  shapeDerivative.resize(order + 1);
  for (auto i = 0; i < 2; ++i) {
    // LobattoShape(shape(i), shapeDerivative(i), xsi, i);
    LobattoShapes::getShape(shape(i), shapeDerivative(i), xsi, i);
  }
  Types::VectorX<prec> tempshape;
  Types::VectorX<prec> tempshapeDerivative;
  this->getH1ShapesInternal(order, tempshape, tempshapeDerivative,
                            xsi);
  for (auto i = 0; i < tempshape.rows(); ++i) {
    shape(i + 2) = tempshape(i);
    shapeDerivative(i + 2) = tempshapeDerivative(i);
  }
}

void LinearEdgeData::getH1ShapesInternal(indexType order,
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

auto LinearEdgeData::getH1Shapes(indexType order,
                                 IntegrationPoint &integration_point)
    -> H1Shapes {
  H1Shapes shapes(order + 1, 1);

  for (auto i = 0; i < 2; ++i) {
    auto sh = LobattoShapes::getShape(integration_point.xi, i);
    shapes.shapes(i) = sh.shapeValue;
    shapes.shapeDeriv(0, i) = sh.shapeDerivative;
  }
  // if (order > 1) {
  auto addShapes =
      this->getH1ShapesInternal(order, integration_point);
  shapes.shapes.tail(order - 1) = addShapes.shapes;
  // shapes.shapes.block(0, 2, 1, order - 1) =
  //     addShapes.shapes.block(0, 0, 1, order - 1);
  shapes.shapeDeriv.block(0, 2, 1, order - 1) =
      addShapes.shapeDeriv.block(0, 0, 1, order - 1);
  //}

  return shapes;
}

auto LinearEdgeData::getH1ShapesInternal(indexType order,
                                         IntegrationPoint &integration_point)
    -> H1Shapes {

  H1Shapes shapes(order - 1, 1);
  for (auto i = 2; i < order + 1; ++i) {
    auto sh = LobattoShapes::getShape(integration_point.xi, i);
    shapes.shapes(i - 2) = sh.shapeValue;
    shapes.shapeDeriv(0, i - 2) = sh.shapeDerivative;
  }

  return shapes;
}

auto LinearEdgeData::getJacobian(prec xi) -> prec {

  prec jac;

  Types::Vector3<prec> dirVec;
  dirVec = m_verts_pointers[1]->getCoordinates();
  dirVec -= m_verts_pointers[0]->getCoordinates();

  jac = dirVec.norm() / prec(2);

  return jac;
}

auto LinearEdgeData::getJacobian(IntegrationPoint &IntegrationPt)
    -> prec {

  prec jacob = 0;
  Types::Vector3<prec> dirVec;
  dirVec = m_verts_pointers[1]->getCoordinates();
  dirVec -= m_verts_pointers[0]->getCoordinates();

  jacob = dirVec.norm() / static_cast<prec>(2);

  return jacob;
}

void LinearEdgeData::setLoad(LoadList &loadlist, indexType meshid,
                             ShapeFunctionTypes shapeType, indexType shapeOrder,
                             Types::VectorX<prec> &Loads, indexType propNumber,
                             Types::VectorX<prec> &direction, bool local,
                             bool add) {
  if (shapeType == ShapeFunctionTypes::HDiv ||
      shapeType == ShapeFunctionTypes::HCurl) {

    std::vector<GenericNodes *> tempNodes;
    this->getNodes(tempNodes, meshid);

    auto tempDofs = tempNodes[0]->getDegreesOfFreedom();
    for (auto i = 0; i < 3; ++i) {
      loadlist.setLoad(propNumber, tempDofs[i]->getId(),
                                      Loads(i), false);
    }

  } else {
    auto GP = this->getIntegrationPoints(-1);
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
    this->getNodes(tempNodes, meshid);

    if (shapeOrder > static_cast<indexType>(tempNodes.size() - 1)) {
      shapeOrder = tempNodes.size() - 1;
    }

    if (!add) {
      for (auto &tempNode : tempNodes) {
        auto tempDofs = tempNode->getDegreesOfFreedom();
        for (auto i = 0; i < 3; ++i) {
          loadlist.setLoad(propNumber, tempDofs[i]->getId(),
                                          prec(0), false);
        }
      }
    }

    for (auto i : GP) {
      auto jacobi = this->getJacobian(i);
      auto h1Shapes = this->getH1Shapes(shapeOrder, i);
      prec dA = jacobi * i.weight;

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
        auto tempDofs = tempNodes[nn]->getDegreesOfFreedom();
        for (auto j = 0; j < 3; ++j) {
          prec loadValue = localLoad(j);
          loadValue *= h1Shapes.shapes(nn);
          loadValue *= dA;
          loadlist.setLoad(propNumber, tempDofs[j]->getId(),
                                          loadValue, true);
        }
      }
    }
  }
}

void LinearEdgeData::setPrescribedSolution(
    LoadList &loadlist, indexType meshid, ShapeFunctionTypes shapeType,
    indexType shapeOrder, Types::VectorX<prec> &Solution, indexType propNumber,
    Types::VectorX<prec> &direction, bool local, bool add) {
  if (shapeType == ShapeFunctionTypes::H1) {
    std::vector<GenericNodes *> tempNodes =
        this->getH1Nodes(meshid, shapeOrder);
    for (indexType i = 0; i < static_cast<indexType>(tempNodes.size()); i++) {
      for (indexType j = 0; j < 3; ++j) {
        if (abs(Solution(j)) > std::numeric_limits<prec>::epsilon() * 1000) {
          auto &dof = tempNodes[i]->getDegreeOfFreedom(j);
          dof.setStatus(dofStatus::inactive);
          if (j < 2)
            loadlist.setLoad(
                propNumber, dof.getId(), Solution(j), add);
        }
      }
    }
  } 
}


auto LinearEdgeData::getDirectionVector()
    -> Types::Vector3<prec> {
  Types::Vector3<prec> dir;
  dir = m_verts_pointers[1]->getCoordinates() -
        m_verts_pointers[0]->getCoordinates();
  return dir;
}

auto LinearEdgeData::getA1Vector(IntegrationPoint &integration_point)
    -> Types::Vector3<prec> {

  Types::Vector3<prec> dir;
  dir = m_verts_pointers[1]->getCoordinates() -
        m_verts_pointers[0]->getCoordinates();
  dir /= dir.norm();
  return dir;
}

void LinearEdgeData::setBoundaryCondition(indexType meshId, indexType order,
                                          ShapeFunctionTypes shapeType,
                                          Types::Vector3<indexType> &dofs,
                                          bool set) {

  if (shapeType == ShapeFunctionTypes::H1) {
    std::vector<GenericNodes *> tnodes;
    this->getNodes(tnodes, meshId);
    indexType nnodes = order + 1;
    if (static_cast<indexType>(tnodes.size()) < nnodes)
      nnodes = tnodes.size();
    if (set) {
      for (auto i = 0; i < nnodes; ++i) {
        GenericNodes *tnode;
        tnode = tnodes[i];
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            tnode->setBoundaryCondition(j);
          } else {
            tnode->unsetBoundaryCondition(j);
          }
        }
      }
    } else {
      for (auto i = 0; i < nnodes; ++i) {
        GenericNodes *tnode;
        tnode = tnodes[i];
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            tnode->setBoundaryCondition(j);
          }
        }
      }
    }

  } else if (shapeType == ShapeFunctionTypes::H0 ||
             shapeType == ShapeFunctionTypes::HCurl ||
             shapeType == ShapeFunctionTypes::HCurl) {
    std::vector<GenericNodes *> tnodes;
    tnodes = this->getNodesOfSet(meshId);
    indexType nnodes = order - 1;
    if (static_cast<indexType>(tnodes.size()) < nnodes)
      nnodes = tnodes.size();
    if (set) {
      for (auto i = 0; i < nnodes; ++i) {
        GenericNodes *tnode;
        tnode = tnodes[i];
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            tnode->setBoundaryCondition(j);
          } else {
            tnode->unsetBoundaryCondition(j);
          }
        }
      }
    } else {
      for (auto i = 0; i < nnodes; ++i) {
        GenericNodes *tnode;
        tnode = tnodes[i];
        for (auto j = 0; j < 3; ++j) {
          if (dofs(j) != 0) {
            tnode->setBoundaryCondition(j);
          }
        }
      }
    }
  }
}

void LinearEdgeData::getNodesInternal(std::vector<GenericNodes *> &nodeVector,
                                      indexType meshId) {
  nodeVector = this->getNodesOfSet(meshId);
}

const GeometryTypes LinearEdgeData::type = GeometryTypes::LinearEdge;

void LinearEdgeData::setAllNodeBoundaryConditionMeshId(
    indexType meshId, indexType dof) {
  GeometryBaseData::setAllNodeBoundaryConditionMeshId(meshId, dof);

  for (auto &i : this->m_verts_pointers) {
    i->setAllNodeBoundaryConditionMeshId(meshId, dof);
  }
}

void LinearEdgeData::flip() {

  std::reverse(this->m_verts.begin(), this->m_verts.end());
  std::reverse(this->m_verts_pointers.begin(), this->m_verts_pointers.end());

}



} // namespace HierAMuS::Geometry
