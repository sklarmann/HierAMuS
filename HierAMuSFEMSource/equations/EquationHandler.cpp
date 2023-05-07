// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include <equations/EquationHandler.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>

#include <equations/DegreeOfFreedom.h>
#include <equations/DofStatus.h>
#include <equations/GenericNodes.h>
#include <equations/NodeSet.h>
#include <pointercollection/pointercollection.h>

#include <solver/GenericSolutionState.h>

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

  indexType sets = 0;
  if (setCollectionIndex < this->m_nodesetCompressedIndex.size() - 1) {
    sets = this->m_nodesetCompressedIndex[setCollectionIndex + 1] -
           this->m_nodesetCompressedIndex[setCollectionIndex];
  }
  return sets;
}

auto EquationHandler::requestNodeSetSetup(indexType setCollectionIndex,
                                          indexType meshID,
                                          bool &alreadyInitialized)
    -> indexType {

  if (setCollectionIndex >= this->m_nodesetCompressedIndex.size() - 1)
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

NodeSet *EquationHandler::getSet(indexType setCollectionIndex,
                                 indexType setNumber) {
  // TODO Adding an exception when trying to access non present node set. //

  if (static_cast<std::size_t>(setCollectionIndex) <
      this->m_nodesetCompressedIndex.size()) {
    indexType pos =
        this->m_nodesetCompressedIndex[setCollectionIndex] + setNumber;
    return &this->m_nodeSetList[pos];
  } else {
    return nullptr;
  }
}

auto EquationHandler::getSetMeshId(indexType setCollectionIndex,
                                   indexType meshId) -> NodeSet * {
  // TODO Adding an exception when trying to access non present node set. //

  indexType sets = this->m_nodesetCompressedIndex[setCollectionIndex + 1] -
                   this->m_nodesetCompressedIndex[setCollectionIndex];
  bool search = true;
  indexType pos = 0;
  NodeSet *tempSet = nullptr;
  indexType sindex = this->m_nodesetCompressedIndex[setCollectionIndex];
  while (pos < sets && search) {
    tempSet = &this->m_nodeSetList[sindex + pos];
    if (tempSet->getMeshId() == meshId) {
      search = false;
    }
    ++pos;
  }

  return tempSet;
}

void EquationHandler::getSets(std::vector<NodeSet *> &nodeSets,
                              indexType setCollectionIndex, indexType setSize) {
  auto colIndex = static_cast<std::size_t>(setCollectionIndex);
  if (colIndex < this->m_nodesetCompressedIndex.size() - 1) {
    indexType setLength =
        this->m_nodesetCompressedIndex[setCollectionIndex + 1] -
        this->m_nodesetCompressedIndex[setCollectionIndex];
    indexType startIndex = this->m_nodesetCompressedIndex[setCollectionIndex];
    nodeSets.clear();
    if (setLength > setSize) {
      // TODO Add exception different sizes.
    }
    for (auto i = 0; i < setLength; ++i) {
      nodeSets.push_back(&this->m_nodeSetList[startIndex + i]);
    }
  } else {
    // TODO Adding an exception when trying to access non present node set
    // collection. //
  }
}

auto EquationHandler::getSets(indexType setCollectionIndex, indexType setSize)
    -> std::vector<NodeSet *> {
  const auto colIndex = static_cast<std::size_t>(setCollectionIndex);
  std::vector<NodeSet *> nodeSets;
  if (colIndex < this->m_nodesetCompressedIndex.size() - 1) {
    const auto setLength =
        this->m_nodesetCompressedIndex[setCollectionIndex + 1] -
        this->m_nodesetCompressedIndex[setCollectionIndex];
    const auto startIndex = this->m_nodesetCompressedIndex[setCollectionIndex];

    if (setLength > setSize) {
      // TODO Add exception different sizes.
    }
    for (auto i = 0; i < setLength; ++i) {
      nodeSets.push_back(&this->m_nodeSetList[startIndex + i]);
    }
  } else {
    // TODO Adding an exception when trying to access non present node set
    // collection. //
  }
  return nodeSets;
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
      this->m_nodes.emplace_back(m_totalIds);        // Add the node to the list
      this->m_nodes.back().setId(this->m_currNodes); // Assign the node Id
      this->m_nodes.back().setNodeType(
          temp->getNodeSetType()); // Assign the node type
      // this->m_nodes.back().requestEquationIds(
      //     *this->pointers); // Request the Specific amount of degrees of
      //     freedom
      ++this->m_currNodes;
    }
    this->m_nodesIndex.push_back(static_cast<indexType>(this->m_nodes.size()));

    ++this->m_currNodeSets;
  }
}

void EquationHandler::print(PointerCollection &pointers) {

  auto &Logger = pointers.getSPDLogger();

  Logger.info("\n{:-<100}\n"
    "Information on EquationHandler:\n"
         "   Total Nodesets:                {:>12}\n"
         "   Total Nodes:                   {:>12}\n"
         "   Total Degrees of Freedom:      {:>12}\n"
         "   Active Degrees of Freedom:     {:>12}\n"
         "   Inactive Degrees of Freedom:   {:>12}\n"
         "\n{:-<100}",
         "",m_nodeSetList.size(),m_nodes.size(),m_totalIds,m_activeIds,m_inActiveIds,"");

}

auto EquationHandler::getNode(const NodeSet &nodeSet, indexType node)
    -> GenericNodes * {
  const auto nodesInSet = nodeSet.getNumberOfNodes();
  if (node >= nodesInSet) {
    // TODO throw exception
    std::cout << "Node not in set" << std::endl;
  }
  const auto storageId = nodeSet.getNodeStorageStartId();
  const auto compressedIndex = this->m_nodesIndex[storageId];
  return &this->m_nodes[compressedIndex + node];
}

void EquationHandler::getNodes(std::vector<GenericNodes *> &nodes,
                               indexType setCollectionIndex,
                               indexType setNumber) {
  if (setCollectionIndex >=
      static_cast<indexType>(this->m_nodesetCompressedIndex.size())) {
    // TODO throw exception
    std::cout << "Node set not present" << std::endl;
  }
  const auto sets = this->m_nodesetCompressedIndex[setCollectionIndex + 1] -
                    this->m_nodesetCompressedIndex[setCollectionIndex];
  if (setNumber >= sets) {
    // TODO throw exception
    std::cout << "Set not in setgroup" << std::endl;
  }
  indexType pos =
      this->m_nodesetCompressedIndex[setCollectionIndex] + setNumber;

  auto tempSet = &this->m_nodeSetList[pos];
  const auto numberOfNodes = tempSet->getNumberOfNodes();
  nodes.resize(numberOfNodes);

  for (auto i = 0; i < numberOfNodes; ++i) {
    nodes[i] = this->getNode(*tempSet, i);
  }
}

auto EquationHandler::getNodes(indexType setCollectionIndex,
                               indexType setNumber)
    -> std::vector<GenericNodes *> {
  if (setCollectionIndex >=
      static_cast<indexType>(this->m_nodesetCompressedIndex.size())) {
    // TODO throw exception
    std::cout << "Node set not present" << std::endl;
  }
  const auto sets = this->m_nodesetCompressedIndex[setCollectionIndex + 1] -
                    this->m_nodesetCompressedIndex[setCollectionIndex];
  if (setNumber >= sets) {
    // TODO throw exception
    std::cout << "Set not in setgroup" << std::endl;
  }
  indexType pos =
      this->m_nodesetCompressedIndex[setCollectionIndex] + setNumber;
  NodeSet *tempSet;
  tempSet = &this->m_nodeSetList[pos];
  const auto numberOfNodes = tempSet->getNumberOfNodes();
  std::vector<GenericNodes *> nodes;
  nodes.resize(numberOfNodes);

  for (auto i = 0; i < numberOfNodes; ++i) {
    nodes[i] = this->getNode(*tempSet, i);
  }
  return nodes;
}

void EquationHandler::getNodes(std::vector<GenericNodes *> &nodesOut,
                               const NodeSet &Nodeset) {
  indexType NodeStorageId = Nodeset.getNodeStorageStartId();
  indexType numOfNodes = Nodeset.getNumberOfNodes();
  nodesOut.clear();

  for (auto i = 0; i < numOfNodes; ++i) {
    nodesOut.push_back(&this->m_nodes[this->m_nodesIndex[NodeStorageId] + i]);
  }
}

auto EquationHandler::getNodes(const NodeSet &Nodeset)
    -> std::vector<GenericNodes *> {
  const auto NodeStorageId = Nodeset.getNodeStorageStartId();
  const auto numOfNodes = Nodeset.getNumberOfNodes();
  std::vector<GenericNodes *> nodesOut;
  nodesOut.reserve(numOfNodes);

  for (auto i = 0; i < numOfNodes; ++i) {
    nodesOut.push_back(&this->m_nodes[this->m_nodesIndex[NodeStorageId] + i]);
  }
  return nodesOut;
}

void EquationHandler::initSolutionState(PointerCollection &pointers) {
  pointers.getSolutionState()->setInitialValues(this->m_totalIds,
                                                this->m_activeIds);
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

} /* namespace HierAMuS */
