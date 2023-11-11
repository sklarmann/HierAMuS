// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "NodeSetManager.h"

#include "EquationHandler.h"
#include "GenericNodes.h"

#include <ostream>

#include <iomanip>
#include <sstream>

namespace HierAMuS {
NodeSetManager::NodeSetManager() : m_nodeSetStorageId(-1), m_numberOfNodeSets(0), m_eqHandler(NULL) {}
NodeSetManager::NodeSetManager(EquationHandler *eqHandler)
    : m_nodeSetStorageId(-1), m_numberOfNodeSets(0), m_eqHandler(eqHandler) {}

NodeSetManager::~NodeSetManager() {}

void NodeSetManager::addNewNodeSet(indexType meshId, indexType numberOfNodes,
                                   NodeTypes type) {

  bool init;
  m_nodeSetStorageId == -1 ? init = false : init = true;
  m_nodeSetStorageId =
      m_eqHandler->requestNodeSetSetup(m_nodeSetStorageId, meshId, init);
  m_numberOfNodeSets = m_eqHandler->getNumberOfNodeSets(m_nodeSetStorageId);
  auto sets = m_eqHandler->getNodeSetList(m_nodeSetStorageId);
  auto &set = sets.getSetMeshId(meshId);
  set.setNumberOfNodes(numberOfNodes);
  set.setType(type);
}

auto NodeSetManager::getNodeSetList() -> NodeSetList {
  return m_eqHandler->getNodeSetList(m_nodeSetStorageId);
}

auto NodeSetManager::getNodeSetNodeList(indexType meshId) -> NodeSetNodeList {
  auto sets = m_eqHandler->getNodeSetList(m_nodeSetStorageId);
  if (sets.hasSetMeshId(meshId)) {
    return m_eqHandler->getNodeSetNodeList(sets.getSetMeshId(meshId));
  }
  NodeSet tempset(NodeTypes::undef, meshId, 0);
  return m_eqHandler->getNodeSetNodeList(tempset);
}

auto NodeSetManager::getNodeVector(indexType meshId)
    -> std::vector<GenericNodes *> {
  
  return std::vector<GenericNodes *>();
}

auto NodeSetManager::getAllDegreesOfFreedom()
    -> std::vector<DegreeOfFreedom *> {
  auto sets = m_eqHandler->getNodeSetList(m_nodeSetStorageId);
  std::vector<DegreeOfFreedom *> Dofs;
  for (auto &set : sets) {
    auto nodes = m_eqHandler->getNodeSetNodeList(set);
    auto temp = nodes.getDegreesOfFreedom();
    Dofs.insert(Dofs.end(), temp.begin(), temp.end());
  }
  return Dofs;
}

void NodeSetManager::print(spdlog::logger &Log) {

  Log.debug("   Number of nodesets:  {:>12}",
               m_numberOfNodeSets);
  if (m_numberOfNodeSets > 0) {
    auto nodeSets = this->getNodeSetList();
    for (auto &set : nodeSets) {
      Log.debug(set);
      auto nodes = m_eqHandler->getNodeSetNodeList(set);
      for (auto &node : nodes) {
        Log.debug(node);
      }
    }
  }
}

} /* namespace HierAMuS */
