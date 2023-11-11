// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <geometry/GeometryBaseData.h>

#include "geometry/GeometryData.h"

#include "EquationHandler.h"

#include "shapefunctions/IntegrationsPoints/helperClasses/IntegrationPointsManagement.h"

namespace HierAMuS::Geometry {

GeometryBaseData::GeometryBaseData() {}

GeometryBaseData::~GeometryBaseData() = default;

auto GeometryBaseData::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::GeometryBaseData::type;
};
auto GeometryBaseData::getGroupType() -> const GeometryTypes & {
  return HierAMuS::Geometry::GeometryBaseData::type;
};

void GeometryBaseData::setNodeSet(indexType meshID, indexType numberOfNodes,
                                  NodeTypes type) {

  m_NodeSetManager.addNewNodeSet(meshID, numberOfNodes, type);
}


void GeometryBaseData::printEqInfo(spdlog::logger &Logger) {

  Logger.debug("   Geometry element's id: {:>12}", this->id);

  m_NodeSetManager.print(Logger);
}

void GeometryBaseData::setAllNodeBoundaryConditionMeshId(indexType meshId,
                                                         indexType dof) {
  auto nodeList = m_NodeSetManager.getNodeSetNodeList(meshId);
  for (auto &node : nodeList) {
    node.setBoundaryCondition(dof);
  }
}

auto GeometryBaseData::getNodesOfSet(indexType meshId)
    -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodeVector;

  auto setList = m_NodeSetManager.getNodeSetNodeList(meshId);
  return setList.getNodes();
}

auto GeometryBaseData::getNodeSetNodeListMeshId(indexType meshId)
    -> NodeSetNodeList {
  return m_NodeSetManager.getNodeSetNodeList(meshId);
}

auto GeometryBaseData::getNodeSetList() -> NodeSetList {
  return m_NodeSetManager.getNodeSetList();
}

void GeometryBaseData::getAllEquationsIds(
    std::vector<DegreeOfFreedom *> &Dofs) {
  Dofs = m_NodeSetManager.getAllDegreesOfFreedom();
}

void GeometryBaseData::getNodeEquationIds(std::vector<DegreeOfFreedom *> &Dofs,
                                          indexType meshId,
                                          indexType nodeNumber) {
  auto nodes = m_NodeSetManager.getNodeSetNodeList(meshId);
  if (nodeNumber < nodes.getNumberOfNodes()) {
    Dofs = nodes[nodeNumber].getDegreesOfFreedom();
  }
}

void GeometryBaseData::getNodesInternal(std::vector<GenericNodes *> &nodeVector,
                                        indexType meshId) {}

auto GeometryBaseData::getNodeListMap()
    -> std::map<indexType, NodeSetNodeList> {

  std::map<indexType, NodeSetNodeList> NodeList;
  auto sets = m_NodeSetManager.getNodeSetList();
  for (auto &it : sets) {
    auto meshId = it.getMeshId();
    NodeList.emplace(meshId, m_NodeSetManager.getNodeSetNodeList(meshId));
  }
  return NodeList;
}

auto GeometryBaseData::getIntegrationPoints(indexType elementId)
    -> IntegrationPoints {
  return IntegrationPointsManagement::getIntegrationsPoints(-1);
}

const HierAMuS::Geometry::GeometryTypes GeometryBaseData::type =
    HierAMuS::Geometry::GeometryTypes::Generic;

} // namespace HierAMuS::Geometry
