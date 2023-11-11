// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include "geometry/GeometryData.h"
#include "geometry/GeometryTypes.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "shapefunctions/LobattoShapes.h"

#include "NodeSet.h"
#include "geometry/Edges/LinearEdgeRuntime.h"
#include "geometry/VertexData.h"

#include "plot/vtkplotClass.h"

#include "DegreeOfFreedom.h"
#include "GenericNodes.h"

#include "LoadList.h"

#include <shapefunctions/LagrangeShape.h>
#include <stdexcept>
#include <types/MatrixTypes.h>
#include <vtkCellType.h>

#include <iomanip>
#include <sstream>

#include "shapefunctions/IntegrationsPoints/helperClasses/IntegrationPointsManagement.h"

namespace HierAMuS::Geometry {
using std::vector;

LinearEdgeRuntime::LinearEdgeRuntime(GeometryData &geoData,
                                     LinearEdgeData &data_element)
    : EdgesRuntimeDataInterface(geoData, data_element) {}

LinearEdgeRuntime::~LinearEdgeRuntime() = default;

auto LinearEdgeRuntime::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::LinearEdgeRuntime::type;
}

void LinearEdgeRuntime::getNodes(std::vector<GenericNodes *> &nodeVector,
                                 indexType meshId) {
  nodeVector.clear();
  for (auto &i : m_Vertices) {
    auto tempNodes = i->getNodesOfSet(meshId);
    nodeVector.insert(nodeVector.end(), tempNodes.begin(), tempNodes.end());
  }
  auto tempNodes = this->getNodesOfSet(meshId);
  nodeVector.insert(nodeVector.end(), tempNodes.begin(), tempNodes.end());
}

auto LinearEdgeRuntime::getH1NodesList(indexType meshID, indexType order)
    -> MeshIdNodeList {
  MeshIdNodeList tempList(meshID);
  tempList.reserve(3);
  for (auto &i : m_Vertices) {
    auto ll = i->getNodeSetNodeListMeshId(meshID);
    tempList.add(ll);
  }
  if (order > 1) {
    auto ll = this->getNodeSetNodeListMeshId(meshID);
    tempList.add(ll);
  }
  return tempList;
}

void LinearEdgeRuntime::print(spdlog::logger &Log) {
  this->m_Data_Element->print(Log);
}

auto LinearEdgeRuntime::getEdgeOrientation(indexType startVertex,
                                           indexType endVertex) -> prec {
  prec val;
  startVertex == m_Vertices[0]->getId() ? val = prec(1.0) : val = prec(-1.0);
  return val;
}

auto LinearEdgeRuntime::getCoordinates(prec xi) -> Types::Vector3<prec> {
  prec v1 = prec(1) - xi;
  v1 /= prec(2);
  prec v2 = prec(1) + xi;
  v2 /= prec(2);

  Types::Vector3<prec> c1;
  Types::Vector3<prec> c2;
  Types::Vector3<prec> c3;

  c1 = m_Vertices[0]->getCoordinates();
  c2 = m_Vertices[1]->getCoordinates();

  c3 = v1 * c1 + v2 * c2;

  return c3;
}

auto LinearEdgeRuntime::getIntegrationPoints(indexType elementId)
    -> IntegrationPoints {

  IntegrationPoints intpoints =
      IntegrationPointsManagement::getIntegrationsPoints(elementId);
  intpoints.setType(IntegrationType::Gauss1D);
  return intpoints;
}

void LinearEdgeRuntime::setH1Shapes(indexType meshId, indexType order,
                                    NodeTypes type) {

  for (auto i = 0; i < 2; ++i) {
    m_Vertices[i]->setNodeSet(meshId, 1, type);
  }
  this->setH1ShapesInternal(meshId, order, type);
}

void LinearEdgeRuntime::setH1ShapesInternal(indexType meshId, indexType order,
                                            NodeTypes type) {
  if (order > 1) {
    this->setNodeSet(meshId, order - 1, type);
  }
}

void LinearEdgeRuntime::getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs,
                                  indexType meshID, indexType order) {

  auto ll = this->getH1NodesList(meshID, order);
  Dofs = ll.getDegreesOfFreedom();
}

void LinearEdgeRuntime::getH1DofsInternal(std::vector<DegreeOfFreedom *> &Dofs,
                                          indexType meshID, indexType order) {

  if (order > 1) {
    if (m_Nodes.find(meshID) != m_Nodes.end()) {
      m_Nodes[meshID].addDofsToVector(Dofs);
    }
  }
}

auto LinearEdgeRuntime::getH1Nodes(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodes;
  for (auto &i : this->m_Vertices) {
    auto tempNodes = i->getNodesOfSet(meshID);
    nodes.insert(nodes.end(), tempNodes.begin(), tempNodes.end());
  }
  if (order > 1) {
    auto tempnodes = this->getH1NodesInternal(meshID, order);
    nodes.insert(nodes.end(), tempnodes.begin(), tempnodes.end());
  }
  return nodes;
}

auto LinearEdgeRuntime::getH1NodesInternal(indexType meshID, indexType order)
    -> std::vector<GenericNodes *> {

  std::vector<GenericNodes *> nodes;
  if (m_Nodes.find(meshID) != m_Nodes.end()) {
    nodes.reserve(m_Nodes[meshID].getNumberOfNodes());
    for (auto &it : m_Nodes[meshID]) {
      nodes.push_back(&it);
    }
  }

  return nodes;
}

auto LinearEdgeRuntime::getH1Shapes(indexType order,
                                    IntegrationPoint &integration_point)
    -> H1Shapes {
  H1Shapes shapes(order + 1, 1);

  for (auto i = 0; i < 2; ++i) {
    auto sh = LobattoShapes::getShape(integration_point.xi, i);
    shapes.shapes(i) = sh.shapeValue;
    shapes.shapeDeriv(0, i) = sh.shapeDerivative;
  }
  // if (order > 1) {
  auto addShapes = this->getH1ShapesInternal(order, integration_point);
  shapes.shapes.tail(order - 1) = addShapes.shapes;
  // shapes.shapes.block(0, 2, 1, order - 1) =
  //     addShapes.shapes.block(0, 0, 1, order - 1);
  shapes.shapeDeriv.block(0, 2, 1, order - 1) =
      addShapes.shapeDeriv.block(0, 0, 1, order - 1);
  //}

  return shapes;
}

auto LinearEdgeRuntime::getH1ShapesInternal(indexType order,
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

void LinearEdgeRuntime::setLoad(
    LoadList &loadlist, indexType meshid, ShapeFunctionTypes shapeType,
    indexType shapeOrder, Types::VectorX<prec> &Loads, indexType propNumber,
    Types::VectorX<prec> &direction, bool local, bool add) {
  if (shapeType == ShapeFunctionTypes::HDiv ||
      shapeType == ShapeFunctionTypes::HCurl) {

    std::vector<GenericNodes *> tempNodes;
    this->getNodes(tempNodes, meshid);

    auto tempDofs = tempNodes[0]->getDegreesOfFreedom();
    for (auto i = 0; i < 3; ++i) {
      loadlist.setLoad(propNumber, tempDofs[i]->getId(), Loads(i), false);
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
          loadlist.setLoad(propNumber, tempDofs[i]->getId(), prec(0), false);
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
          loadlist.setLoad(propNumber, tempDofs[j]->getId(), loadValue, true);
        }
      }
    }
  }
}

void LinearEdgeRuntime::setPrescribedSolution(
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
            loadlist.setLoad(propNumber, dof.getId(), Solution(j), add);
        }
      }
    }
  }
}

void LinearEdgeRuntime::geometryToParaview(vtkPlotInterface &paraviewAdapter,
                                           indexType mainMesh,
                                           indexType subMesh) {
  std::vector<indexType> points(2);
  points = {m_Vertices[0]->getId(), m_Vertices[1]->getId()};
  for (auto &i : this->m_Vertices) {
    i->geometryToParaview(paraviewAdapter, mainMesh, subMesh);
  }
  paraviewAdapter.addCell(mainMesh, subMesh, m_Data_Element->getId(), 1, points,
                          2, VTK_LINE);
}

void LinearEdgeRuntime::computeWeightsParaview(
    vtkPlotInterface &paraviewAdapter, indexType mainMesh, indexType subMesh) {
  auto GP = this->getIntegrationPoints(-1);
  GP.setOrder(2);
  for (auto i : GP) {
    auto jaco = this->getJacobian(i);
    auto shapes = this->getH1Shapes(1, i);
    prec dA = jaco * i.weight;
    for (auto j = 0; j < 2; ++j) {
      std::vector<prec> val;
      val.push_back(shapes.shapes(j) * dA);
      paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, val,
                                           m_Vertices[j]->getId(), 1,
                                           paraviewNames::weightName());
    }
  }
}

void LinearEdgeRuntime::H1SolutionToParaview(
    vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    Types::VectorX<prec> solution, std::string &name) {

  for (auto i = 0; i < 2; ++i) {
    std::vector<prec> sol(3);
    for (auto j = 0; j < 3; ++j) {
      sol[j] = solution(3 * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, m_Vertices[i]->getId(), sol,
                                 3, name);
  }
}

void LinearEdgeRuntime::H1DataToParaview(vtkPlotInterface &paraviewAdapter,
                                         indexType mainMesh, indexType subMesh,
                                         Types::VectorX<prec> &Data,
                                         indexType numberComponents,
                                         indexType order, std::string &name) {
  for (auto i = 0; i < 2; ++i) {
    std::vector<prec> sol(numberComponents);
    for (auto j = 0; j < numberComponents; ++j) {
      sol[j] = Data(numberComponents * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, m_Vertices[i]->getId(), sol,
                                 numberComponents, name);
  }
}

void LinearEdgeRuntime::projectDataToParaviewVertices(
    vtkPlotInterface &paraviewAdapter, indexType mainMesh, indexType subMesh,
    indexType order, IntegrationPoint &IntegrationPt,
    Types::VectorX<prec> &data, indexType numberComponents, std::string name) {
  auto shapes = this->getH1Shapes(1, IntegrationPt);
  std::vector<prec> vals(numberComponents);

  auto jaco = this->getJacobian(IntegrationPt);
  auto dA = jaco * IntegrationPt.weight;
  for (auto i = 0; i < 2; ++i) {
    for (auto j = 0; j < numberComponents; ++j) {
      vals[j] = data(j) * shapes.shapes(i) * dA;
    }
    paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, vals,
                                         m_Vertices[i]->getId(),
                                         numberComponents, name);
  }
}

auto LinearEdgeRuntime::getDirectionVector() -> Types::Vector3<prec> {
  Types::Vector3<prec> dir;
  dir = m_Vertices[1]->getCoordinates() - m_Vertices[0]->getCoordinates();
  return dir;
}

auto LinearEdgeRuntime::getA1Vector(IntegrationPoint &integration_point)
    -> Types::Vector3<prec> {

  Types::Vector3<prec> dir;
  dir = m_Vertices[1]->getCoordinates() - m_Vertices[0]->getCoordinates();
  dir /= dir.norm();
  return dir;
}

void LinearEdgeRuntime::setBoundaryCondition(indexType meshId, indexType order,
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

void LinearEdgeRuntime::getNodesInternal(
    std::vector<GenericNodes *> &nodeVector, indexType meshId) {
  nodeVector = this->getNodesOfSet(meshId);
}

const GeometryTypes LinearEdgeRuntime::type = GeometryTypes::LinearEdge;

void LinearEdgeRuntime::setAllNodeBoundaryConditionMeshId(indexType meshId,
                                                          indexType dof) {
  GeometryBaseRuntime::setAllNodeBoundaryConditionMeshId(meshId, dof);

  for (auto &i : this->m_Vertices) {
    i->setAllNodeBoundaryConditionMeshId(meshId, dof);
  }
}

} // namespace HierAMuS::Geometry
