// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once



#include <ostream>
#include <vector>
#include "NodeSetList.h"
#include "NodeSetNodeList.h"
#include "MeshIdNodeList.h"
#include "GenericNodes.h"
#include "NodeSetManager.h"

#include <vector>

namespace HierAMuS {
    /**
 * @class EquationHandler
 * @brief Groups nodeSets and dofHandler for further management of the equation
 * system.
 */

class EquationHandler {
public:
  explicit EquationHandler();
  ~EquationHandler();
  auto requestNodeSetSetup(indexType setCollectionIndex, indexType meshID,
                           bool &alreadyInitialized) -> indexType;
  /**
   * @brief Request a node set collection.
   * @param numberOfSets Specifies the number of node sets.
   * @return Returns the storage number of the node set collection to
   * identify the set collection later on.
   */
  // auto requestNodeSets(indexType numberOfSets) -> indexType;
  auto getNumberOfNodeSets(indexType setCollectionIndex) -> indexType;


  void update();
  void updateEquations(); // TODO

  
  auto getNumberOfTotalEquations() -> indexType;
  auto getNumberOfActiveEquations() -> indexType;
  auto getNumberOfInActiveEquations() -> indexType;

  auto getDegreeOfFreedom(indexType dofId) -> DegreeOfFreedom &;
  auto getNumberOfNodes() -> indexType;
  auto getNode(indexType globalNodeNumber) -> GenericNodes &;


  auto getNodeSetNodeList(NodeSet &Set) -> NodeSetNodeList;
  auto getNodeSetList(indexType setCollectionIndx) -> NodeSetList;

  auto getNewNodeSetManager() -> NodeSetManager {
    return NodeSetManager(this);
  };

  template <typename OStream>
  friend OStream &operator<<(OStream &os, const EquationHandler &self) {
    fmt::format_to(std::ostream_iterator<char>(os),
                   "\n{:-<100}\n"
                   "Information on EquationHandler:\n"
                   "   Total Nodesets:                {:>12}\n"
                   "   Total Nodes:                   {:>12}\n"
                   "   Total Degrees of Freedom:      {:>12}\n"
                   "   Active Degrees of Freedom:     {:>12}\n"
                   "   Inactive Degrees of Freedom:   {:>12}\n"
                   "\n{:-<100}",
                   "", self.m_nodeSetList.size(), self.m_nodes.size(),
                   self.m_totalIds, self.m_activeIds, self.m_inActiveIds, "");
    return os;
  }


private:
  std::vector<indexType> m_nodesetCompressedIndex;
  std::vector<NodeSet> m_nodeSetList;
  indexType m_currNodeSets;

  indexType m_currNodes;
  std::vector<indexType> m_nodesIndex;
  std::vector<GenericNodes> m_nodes;
  

  indexType m_activeIds;
  indexType m_inActiveIds;
  indexType m_totalIds;
};

} /* namespace HierAMuS */
