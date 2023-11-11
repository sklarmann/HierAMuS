// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "datatypes.h"


#include "Nodetypes.h"
#include <vector>

#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
namespace HierAMuS {
class DegreeOfFreedom;
class GenericNodes;
 /**
 * @class NodeSetFinal
 * @brief NodeSetFinal groups  a set of nodes and is associated to a geometric
 * element.
 */

class NodeSetNodeList {
public:
  NodeSetNodeList() = default;
  NodeSetNodeList(const NodeSetNodeList &nodelist);
  NodeSetNodeList(NodeSetNodeList &&nodelist) = default;
  NodeSetNodeList(NodeTypes type, indexType meshId, indexType numberOfNodes,
               std::vector<GenericNodes>::iterator &startNode,
               std::vector<GenericNodes>::iterator &endNode)
      : m_type(type), m_meshId(meshId), m_numberOfNodes(numberOfNodes), m_startNode(startNode),
        m_endNode(endNode){};
  NodeSetNodeList(NodeTypes type, indexType meshId, indexType numberOfNodes,
               std::vector<GenericNodes>::iterator &&startNode,
               std::vector<GenericNodes>::iterator &&endNode)
      : m_type(type), m_meshId(meshId), m_numberOfNodes(numberOfNodes), m_startNode(startNode),
        m_endNode(endNode){};
  ~NodeSetNodeList();

  auto operator=(NodeSetNodeList &other) -> NodeSetNodeList&;

  auto getNodes() -> std::vector<GenericNodes*>;
  auto getDegreesOfFreedom() -> std::vector<DegreeOfFreedom *>;
  void addDofsToVector(std::vector<DegreeOfFreedom *> &Dofs);
  auto getNumberOfNodes() -> indexType;

  auto getMeshId() -> indexType;
  
  auto operator[](indexType indx) -> GenericNodes &;

  auto begin() -> std::vector<GenericNodes>::iterator;
  auto end() -> std::vector<GenericNodes>::iterator;

private:
  NodeTypes m_type;          /**< Node type of the set.  */
  indexType m_meshId;        /**< Used for identification of the NodeSetFinals in the
                              preprocessing steps */
  indexType m_numberOfNodes;
  std::vector<GenericNodes>::iterator m_startNode;
  std::vector<GenericNodes>::iterator m_endNode;
};

} /* namespace HierAMuS */
