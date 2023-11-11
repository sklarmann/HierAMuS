// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "EquationHandler.h"


#include "DegreeOfFreedom.h"
#include "DofStatus.h"
#include "GenericNodes.h"
#include "NodeSet.h"

#include "NodeSetNodeList.h"
#include "NodeSetList.h"

#include <iomanip>
#include <iostream>

namespace HierAMuS {

EquationHandler::EquationHandler()
    : m_currNodeSets(0), m_currNodes(0), m_activeIds(0), m_inActiveIds(0),
      m_totalIds(0) {
  this->m_nodesetCompressedIndex.push_back(0);
  this->m_nodesIndex.push_back(0);
}

EquationHandler::~EquationHandler() {}

auto EquationHandler::getNumberOfNodeSets(indexType setCollectionIndex)
    -> indexType {

  if (setCollectionIndex < static_cast<indexType>(this->m_nodesetCompressedIndex.size() - 1)) {
    return this->m_nodesetCompressedIndex[setCollectionIndex + 1] -
           this->m_nodesetCompressedIndex[setCollectionIndex];
  }
  return 0;
}

auto EquationHandler::requestNodeSetSetup(indexType setCollectionIndex,
                                          indexType meshID,
                                          bool &alreadyInitialized)
    -> indexType {

  if (setCollectionIndex >= static_cast<indexType>(this->m_nodesetCompressedIndex.size() - 1))
    alreadyInitialized = false;

  if (!alreadyInitialized) {
    indexType ret;
    ret = this->m_nodesetCompressedIndex.size() - 1;
    // if(this->nodesetCompressedIndex.size()==1){
    //   ret = 0;
    // } else {
    //   ret = this->nodesetCompressedIndex.back();
    // }
    this->m_nodeSetList.emplace_back();
    this->m_nodeSetList.back().setMeshId(meshID);
    this->m_nodesetCompressedIndex.emplace_back(
        static_cast<indexType>(this->m_nodeSetList.size()));
    alreadyInitialized = true;
    return ret;
  } else {
    indexType pos = this->m_nodesetCompressedIndex[setCollectionIndex];
    indexType end = this->m_nodesetCompressedIndex[setCollectionIndex + 1];
    bool search = true;
    while (pos < end && search) {
      if (this->m_nodeSetList[pos].getMeshId() == meshID)
        search = false;
      ++pos;
    }
    if (search) {
      NodeSet toInsert;
      toInsert.setMeshId(meshID);
      auto it = this->m_nodeSetList.begin() + pos;
      this->m_nodeSetList.insert(it, toInsert);

      for (auto cit =
               this->m_nodesetCompressedIndex.begin() + setCollectionIndex + 1;
           cit != this->m_nodesetCompressedIndex.end(); ++cit) {
        ++*cit;
      }
    }
    return setCollectionIndex;
  }
}



/**
 * @brief Updates the nodes correspondig to node sets.
 */

void EquationHandler::update() {
  auto setSize = static_cast<indexType>(this->m_nodeSetList.size());
  while (this->m_currNodeSets < setSize) {
    NodeSet *temp = &this->m_nodeSetList[m_currNodeSets];

    indexType num;
    num = temp->getNumberOfNodes();
    temp->setNodeStorageId(static_cast<indexType>(this->m_nodesIndex.size()) -
                           1); // Set nodeset storage id
    for (auto i = 0; i < num; ++i) {
      this->m_nodes.emplace_back(m_totalIds, this->m_currNodes); // Add the node to the list
      m_totalIds += 3;
      this->m_nodes.back().setNodeType(
          temp->getType()); // Assign the node type
      ++this->m_currNodes;
    }
    this->m_nodesIndex.push_back(static_cast<indexType>(this->m_nodes.size()));

    ++this->m_currNodeSets;
  }
}



void EquationHandler::updateEquations() {

  m_activeIds = 0;
  m_inActiveIds = 0;
  for (auto &i : m_nodes) {
    for (auto j = 0; j < 3; ++j) {
      auto &dof = i.getDegreeOfFreedom(j);
      if (dof.getStatus() == dofStatus::active) {
        dof.setEqId(m_activeIds);
        ++m_activeIds;
      } else {
        dof.setEqId(m_inActiveIds);
        ++m_inActiveIds;
      }
    }
  }
}

auto EquationHandler::getNumberOfTotalEquations() -> indexType {
  return this->m_totalIds;
}
auto EquationHandler::getNumberOfActiveEquations() -> indexType {
  return this->m_activeIds;
}
auto EquationHandler::getNumberOfInActiveEquations() -> indexType {
  return this->m_inActiveIds;
}

auto EquationHandler::getDegreeOfFreedom(indexType dofId) -> DegreeOfFreedom & {
  indexType masterNodeId = dofId / 3;
  indexType localDofId = dofId % 3;
  return m_nodes[masterNodeId].getDegreeOfFreedom(localDofId);
}

auto EquationHandler::getNumberOfNodes() -> indexType { return m_nodes.size(); }

auto EquationHandler::getNode(indexType globalNodeNumber) -> GenericNodes & {
  return m_nodes[globalNodeNumber];
}



auto EquationHandler::getNodeSetNodeList(NodeSet &Set)
    -> NodeSetNodeList {

  auto storageId = Set.getNodeStorageStartId();
  auto real = true;
  if (storageId < 0) {
    //return NodeSetNodeList(Set.getType(), Set.getMeshId(),
    //                       Set.getNumberOfNodes(), this->m_nodes.end(),
    //                       this->m_nodes.end());
    real = false;
    storageId = 0;
  }



  const auto startIndex = this->m_nodesIndex[storageId];
  const auto endIndex = this->m_nodesIndex[storageId+1];

  auto it1 = this->m_nodes.begin() + startIndex;
  auto it2 = this->m_nodes.begin() + endIndex;

  if (!real) {
    it1 = m_nodes.end();
    it2 = m_nodes.end();
  }

  return NodeSetNodeList(Set.getType(), Set.getMeshId(),
                         Set.getNumberOfNodes(), it1, it2);
}

auto EquationHandler::getNodeSetList(indexType setCollectionIndx)
    -> NodeSetList {
  if(setCollectionIndx == -1){
    return NodeSetList(0, this->m_nodeSetList.end(), this->m_nodeSetList.end());
  }
  auto startIndx = this->m_nodesetCompressedIndex[setCollectionIndx];
  auto endIndx = this->m_nodesetCompressedIndex[setCollectionIndx+1];

  auto it1 = this->m_nodeSetList.begin() + startIndx;
  auto it2 = this->m_nodeSetList.begin() + endIndx;

  return NodeSetList(endIndx - startIndx, it1, it2);
}

} /* namespace HierAMuS */
