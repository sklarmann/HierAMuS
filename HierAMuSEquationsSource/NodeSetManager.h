// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <datatypes.h>

#include "NodeSetList.h"
#include "NodeSetNodeList.h"


#include "Nodetypes.h"
#include <vector>

#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
namespace HierAMuS {
class EquationHandler;
/**
 * @class NodeSetManager
 * @brief NodeSetManager groups  a set of nodes and is associated to a geometric
 * element.
 */

class NodeSetManager {
public:
  NodeSetManager();
  NodeSetManager(EquationHandler *eqHandler);
  NodeSetManager(const NodeSetManager &other) = default;
  NodeSetManager(NodeSetManager &&other) = default;
  ~NodeSetManager();
  
  void addNewNodeSet(indexType meshId, indexType numberOfNodes, NodeTypes type);
  auto getNodeSetList() -> NodeSetList;
  auto getNodeSetNodeList(indexType meshId) -> NodeSetNodeList;
  auto getNodeVector(indexType meshId) -> std::vector<GenericNodes *>;
  auto getAllDegreesOfFreedom() -> std::vector<DegreeOfFreedom *>;


  NodeSetManager &operator=(NodeSetManager &other) = default;
  NodeSetManager &operator=(NodeSetManager &&other) = default;

  void print(spdlog::logger &Log);

private:
  EquationHandler *m_eqHandler;

  indexType m_numberOfNodeSets;
  indexType m_nodeSetStorageId;


};

} /* namespace HierAMuS */
