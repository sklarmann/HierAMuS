// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "datatypes.h"

#include "Nodetypes.h"
#include "NodeSetNodeList.h"

#include <vector>

#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
namespace HierAMuS {

class GenericNodes;
/**
 * @class NodeSetFinalList
 * @brief NodeSetFinalList groups  a set of nodes and is associated to a geometric
 * element.
 */

class MeshIdNodeList {
private:
  struct set {
    using nodeiter = std::vector<GenericNodes>::iterator;
    nodeiter m_start, m_end;
    indexType m_numberOfNodes;

    set(nodeiter start, nodeiter end, indexType numberOfNodes);

    auto begin() -> nodeiter;
    auto end() -> nodeiter;
  };

public:
  MeshIdNodeList() = delete;
  MeshIdNodeList(indexType meshId);
  ~MeshIdNodeList();

  void reserve(indexType size);
  void add(NodeSetNodeList list);
  void append(MeshIdNodeList list);
  auto getDegreesOfFreedom() -> std::vector<DegreeOfFreedom *>;
  auto getNodeVector() -> std::vector<GenericNodes *>;
  auto operator[](indexType indx) -> GenericNodes &;




  struct iterator {
    iterator(std::vector<set>::iterator start, std::vector<set>::iterator end);

    iterator & operator++();

    
    auto operator*() -> GenericNodes &;
    auto operator->() -> GenericNodes *;

    friend bool operator==(const iterator &lhs, const iterator &rhs);
    friend bool operator!=(const iterator &lhs, const iterator &rhs);

    private:
    void update();

    std::vector<set>::iterator m_outerIt, m_outerEnd;
    std::vector<GenericNodes>::iterator m_innerIt, m_innerEnd;
  };

  auto begin() -> iterator;
  auto end() -> iterator;

private:
  void add(std::vector<GenericNodes>::iterator startIter,
           std::vector<GenericNodes>::iterator endIter,
           indexType numberOfNodes);
  
  indexType m_meshId;
  indexType m_numberOfNodes;
  std::vector<set> m_data;
};

} /* namespace HierAMuS */
