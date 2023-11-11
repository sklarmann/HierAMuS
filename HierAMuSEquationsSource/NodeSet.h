// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <datatypes.h>



#include "Nodetypes.h"
#include <vector>

#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
namespace HierAMuS {

/**
 * @class NodeSet
 * @brief Nodeset groups  a set of nodes and is associated to a geometric
 * element.
 */

class NodeSet {
public:
  NodeSet();
  NodeSet(NodeTypes type, indexType meshId, indexType numberOfNodes);
  ~NodeSet();
  /**
   * @brief Set the node type of the set.
   * @param type Node type of the set
   */
  void setType(NodeTypes type);
  /**
   * @brief Returns the type of the set.
   * @return Type of Node.
   */
  auto getType() -> NodeTypes { return this->type; };
  /**
  * @brief Sets the mesh id of the nodes.
  * @param indexType The mesh id of the nodes
  */
  void setMeshId(indexType meshId);
  /**
  * @brief Returns the mesh id of the nodes in the set.
  * @return indexTpye Mesh id.
  */
  auto getMeshId() -> indexType;
  /**
  * @brief Sets the number of nodes the set has.
  * @param numberOfNodes the number of nodes.
  */
  void setNumberOfNodes(indexType numberOfNodes);
  /**
  * @brief Return the number of nodes belonging to the set.
  * @return indexType The number of nodes.
  */
  auto getNumberOfNodes() const -> indexType {
    return this->numberOfNodes;
  };

  /**
  * @brief Sets the storage id of the start nodes in the node list. All nodes of the set are stored continuously.
  * @param indexTpye The storage id.
  */
  void setNodeStorageId(indexType storageId) {
    this->nodeStorageId = storageId;
  };
  /**
  * @brief Returns the storage id of the first node belonging to the set.
  * @return indexType The storage id.
  */
  auto getNodeStorageStartId() const -> indexType {
    return this->nodeStorageId;
  };




  template <typename OStream>
  friend OStream &operator<<(OStream &os, const NodeSet &c) {
    std::string fmtString;
    fmt::format_to(std::ostream_iterator<char>(os),
                   "Nodeset meshId:  {:>5},  Number of Nodes:  {:>5}", c.meshId,
                   c.numberOfNodes);
    return os;
  }

private:
  NodeTypes type;          /**< Node type of the set.  */
  indexType meshId;        /**< Used for identification of the nodeSets in the
                              preprocessing steps */
  indexType numberOfNodes; /**< Stores the number of nodes in the set. */
  indexType
      nodeStorageId; /**< Stores the Storage id of the nodes in the set. */
};

} /* namespace HierAMuS */
