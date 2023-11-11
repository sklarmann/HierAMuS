// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "NodeSetNodeList.h"

#include "EquationHandler.h"
#include "GenericNodes.h"


#include <ostream>

#include <iomanip>
#include <sstream>

namespace HierAMuS {

NodeSetNodeList::NodeSetNodeList(const NodeSetNodeList &other)
    : m_type(other.m_type), m_meshId(other.m_meshId),
      m_numberOfNodes(other.m_numberOfNodes), m_startNode(other.m_startNode),
      m_endNode(other.m_endNode) {}

//NodeSetNodeList::NodeSetNodeList(NodeSetNodeList &&other)
//    : m_type(other.m_type), m_meshId(other.m_meshId),
//      m_numberOfNodes(other.m_numberOfNodes), m_startNode(other.m_startNode),
//      m_endNode(other.m_endNode) {}

NodeSetNodeList::~NodeSetNodeList() = default;

auto NodeSetNodeList::operator=(NodeSetNodeList &other) -> NodeSetNodeList & {
  m_type = other.m_type;
  m_meshId = other.m_meshId;
  m_numberOfNodes = other.m_numberOfNodes;
  m_startNode = other.m_startNode;
  m_endNode = other.m_endNode;
  return *this;
}

auto NodeSetNodeList::getNodes() -> std::vector<GenericNodes *> {

  std::vector<GenericNodes *> Nodes;
  Nodes.reserve(this->m_numberOfNodes);
  for (auto it = this->m_startNode; it != this->m_endNode; it++) {
    Nodes.push_back(&(*it));
  }

  return Nodes;
}

auto NodeSetNodeList::getDegreesOfFreedom() -> std::vector<DegreeOfFreedom *>
{
  std::vector<DegreeOfFreedom *> dofs;
  dofs.reserve(3 * this->m_numberOfNodes);
  for (auto it = this->m_startNode; it != this->m_endNode; ++it) {
    auto tdofs = it->getDegreesOfFreedom();
    dofs.insert(dofs.end(), tdofs.begin(), tdofs.end());
  }
  return dofs;
}

void NodeSetNodeList::addDofsToVector(
    std::vector<DegreeOfFreedom *> &Dofs) {
  for (auto it = this->m_startNode; it != this->m_endNode; ++it) {
    it->addDofsToVector(Dofs);
  }
}

auto NodeSetNodeList::getNumberOfNodes() -> indexType {
  return this->m_numberOfNodes;
}

auto NodeSetNodeList::getMeshId() -> indexType { return m_meshId; }

auto NodeSetNodeList::operator[](indexType indx) -> GenericNodes & {

  if (indx < this->m_numberOfNodes) {
    return *(this->m_startNode + indx);
  }
  throw std::runtime_error(
      "Error in NodeSetNodeList::operator[]: Index out of bounds!");
}

auto NodeSetNodeList::begin() -> std::vector<GenericNodes>::iterator {
  return this->m_startNode;
}

auto NodeSetNodeList::end() -> std::vector<GenericNodes>::iterator {
  return this->m_endNode;
}

} /* namespace HierAMuS */
