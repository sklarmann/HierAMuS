// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <datatypes.h>

#include <forwarddeclaration.h>

#include <equations/Nodetypes.h>
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
  ~NodeSet();
  void setType(NodeTypes type);
  void setMeshId(indexType meshId);
  indexType getMeshId();
  void setNumberOfNodes(indexType numberOfNodes);
  [[nodiscard]] auto getNumberOfNodes() const -> indexType {
    return this->numberOfNodes;
  };
  [[nodiscard]] auto getNodeSetType() const -> NodeTypes { return this->type; };

  void getNodes(PointerCollection &pointers,
                std::vector<GenericNodes *> &Nodes);
  auto getNodes(PointerCollection &pointers) -> std::vector<GenericNodes *>;
  auto getDegreesOfFreedom(PointerCollection &pointers)
      -> std::vector<DegreeOfFreedom *>;

  void setNodeStorageId(indexType storageId) {
    this->nodeStorageId = storageId;
  };
  [[nodiscard]] auto getNodeStorageStartId() const -> indexType {
    return this->nodeStorageId;
  };

  void print(PointerCollection &pointers);

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
