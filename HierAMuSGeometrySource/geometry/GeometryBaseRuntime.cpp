// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "geometry/GeometryBaseData.h"
#include <geometry/GeometryBaseRuntime.h>

#include "GenericNodes.h"
#include "EquationHandler.h"

#include "geometry/GeometryData.h"

#include "shapefunctions/IntegrationsPoints/helperClasses/IntegrationPointsManagement.h"


namespace HierAMuS::Geometry {

GeometryBaseRuntime::GeometryBaseRuntime(GeometryBaseData &data_element)
    : m_Data_Element(&data_element) {

  m_Nodes = data_element.getNodeListMap();
}

GeometryBaseRuntime::~GeometryBaseRuntime() = default;

GeometryBaseRuntime::GeometryBaseRuntime(GeometryBaseRuntime &other)
    : m_Nodes(other.m_Nodes), m_Data_Element(other.m_Data_Element) {}

GeometryBaseRuntime::GeometryBaseRuntime(GeometryBaseRuntime &&other)
    : m_Nodes(std::move(other.m_Nodes)), m_Data_Element(other.m_Data_Element) {}

GeometryBaseRuntime &
GeometryBaseRuntime::operator=(GeometryBaseRuntime &other) {
  m_Nodes = other.m_Nodes;
  m_Data_Element = other.m_Data_Element;
  return *this;
}

GeometryBaseRuntime &
GeometryBaseRuntime::operator=(GeometryBaseRuntime &&other) {
  m_Nodes = std::move(other.m_Nodes);
  m_Data_Element = other.m_Data_Element;
  return *this;
}

void GeometryBaseRuntime::set_id(indexType id) { m_Data_Element->set_id(id); }

auto GeometryBaseRuntime::getId() -> indexType {
  return m_Data_Element->getId();
}

auto GeometryBaseRuntime::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::GeometryBaseRuntime::type;
};
auto GeometryBaseRuntime::getGroupType() -> const GeometryTypes & {
  return HierAMuS::Geometry::GeometryBaseRuntime::type;
};

void GeometryBaseRuntime::setNodeSet(indexType meshID, indexType numberOfNodes,
                                     NodeTypes type) {
  m_Data_Element->setNodeSet(meshID, numberOfNodes, type);
}


void GeometryBaseRuntime::printEqInfo(spdlog::logger &Log) {

  m_Data_Element->printEqInfo(Log);
}

void GeometryBaseRuntime::setAllNodeBoundaryConditionMeshId(indexType meshId,
                                                            indexType dof) {
  if (m_Nodes.find(meshId) != m_Nodes.end()) {
    for (auto &node : m_Nodes[meshId]) {
      node.setBoundaryCondition(dof);
    }
  }
}

auto GeometryBaseRuntime::getNodesOfSet(indexType meshId)
    -> std::vector<GenericNodes *> {

  std::vector<GenericNodes *> ret;
  if (m_Nodes.find(meshId) != m_Nodes.end()) {
    ret.reserve(m_Nodes[meshId].getNumberOfNodes());
    for (auto &it : m_Nodes[meshId]) {
      ret.push_back(&it);
    }
  }
  return ret;
}

auto GeometryBaseRuntime::getNodeSetNodeListMeshId(indexType meshId)
    -> NodeSetNodeList {
  if (m_Nodes.find(meshId) != m_Nodes.end()) {
    return m_Nodes[meshId];
  } else {
    std::vector<GenericNodes> temp;
    temp.clear();
    return NodeSetNodeList(NodeTypes::undef, meshId, 0, temp.end(), temp.end());
  }
}

void GeometryBaseRuntime::getAllEquationsIds(
    std::vector<DegreeOfFreedom *> &Dofs) {
  Dofs.clear();
  for (auto &it : m_Nodes) {
    for (auto &node : it.second) {
      auto nodeDofs = node.getDegreesOfFreedom();
      Dofs.insert(Dofs.end(), nodeDofs.begin(), nodeDofs.end());
    }
  }
}

void GeometryBaseRuntime::getNodeEquationIds(
    std::vector<DegreeOfFreedom *> &Dofs, indexType meshId,
    indexType nodeNumber) {
  if (m_Nodes.find(meshId) != m_Nodes.end()) {
    Dofs = m_Nodes[meshId][nodeNumber].getDegreesOfFreedom();
  }
}

auto GeometryBaseRuntime::checkSetWithMeshId(indexType meshId)
    -> bool {
  if (m_Nodes.find(meshId) != m_Nodes.end()) {
    return true;
  }
  return false;
}

void GeometryBaseRuntime::getNodes(std::vector<GenericNodes *> &nodeVector,
                                   indexType meshId) {}

void GeometryBaseRuntime::getNodesInternal(
    std::vector<GenericNodes *> &nodeVector, indexType meshId) {}

auto GeometryBaseRuntime::getIntegrationPoints(indexType elementId) -> IntegrationPoints {
  return IntegrationPointsManagement::getIntegrationsPoints(elementId);
}

void GeometryBaseRuntime::updateNodes(
    std::map<indexType, NodeSetNodeList> &Nodes) {
  m_Nodes = Nodes;
}

void GeometryBaseRuntime::updateNodes(
    std::map<indexType, NodeSetNodeList> &&Nodes) {
  m_Nodes = std::move(Nodes);
}

const HierAMuS::Geometry::GeometryTypes GeometryBaseRuntime::type =
    HierAMuS::Geometry::GeometryTypes::Generic;

} // namespace HierAMuS::Geometry
