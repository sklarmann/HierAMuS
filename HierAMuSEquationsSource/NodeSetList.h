// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <datatypes.h>

#include "NodeSet.h"

#include "Nodetypes.h"
#include <vector>

#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
namespace HierAMuS {

/**
 * @class NodeSetList
 * @brief NodeSetList groups a set of Nodesets and is associated to a geometric
 * element.
 */

class NodeSetList {
public:
  NodeSetList() = delete;
  /**
  * @brief Constructor for the NodeSetList.
  * @param numberOfSets is the number of NodeSets belonging to the geometric element.
  * @param startIterator Iterator to the first NodeSet belonging to the geometric element.
  * @param endIterator Iterator to the last+1 NodeSet belonging to the geometric element.
  */
  NodeSetList(indexType numberOfSets,
              std::vector<NodeSet>::iterator &startIterator,
              std::vector<NodeSet>::iterator &endIterator);
  NodeSetList(indexType numberOfSets,
              std::vector<NodeSet>::iterator &&startIterator,
              std::vector<NodeSet>::iterator &&endIterator);
  ~NodeSetList();

  /**
  * @brief Operator to access the NodeSets by an integer value.
  */
  auto operator[](indexType indx) -> NodeSet &;
  /**
  * @brief Checks if the List contains a NodeSet with a specific mesh id.
  * @param meshId The mesh id to check for.
  * @return bool True if mesh id is present, false otherwise.
  */
  auto hasSetMeshId(indexType meshId) -> bool;
  /**
  * @brief Returns the set with the specific mesh id.
  * @param meshId The mesh id which the NodeSet must have.
  * @return NodeSet A reference to the NodeSet with the mesh id meshId.
  */
  auto getSetMeshId(indexType meshId) -> NodeSet &;

  /**
  * @brief Returns the number of NodeSets in the list.
  */
  auto getNumberOfSets() -> indexType;

  /**
  * @brief Returns the iterator to the first NodeSet of the list.
  */
  auto begin() -> std::vector<NodeSet>::iterator;
  /**
  * @brief Returns the iterator to the last NodeSet of the list.
  */
  auto end() -> std::vector<NodeSet>::iterator;

private:
  indexType numberOfNodeSets;

  std::vector<NodeSet>::iterator startIterator;
  std::vector<NodeSet>::iterator endIterator;
};

} /* namespace HierAMuS */
