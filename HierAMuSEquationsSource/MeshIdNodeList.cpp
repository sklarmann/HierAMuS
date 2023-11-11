// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "MeshIdNodeList.h"

#include "GenericNodes.h"

#include <iomanip>
#include <sstream>

namespace HierAMuS {

MeshIdNodeList::MeshIdNodeList(indexType meshId)
    : m_meshId(meshId), m_numberOfNodes(0) {}

MeshIdNodeList::~MeshIdNodeList() {}

void MeshIdNodeList::reserve(indexType size) { m_data.reserve(size); }

void MeshIdNodeList::add(NodeSetNodeList list) {
  if (list.getMeshId() == m_meshId) {
    m_data.emplace_back(list.begin(), list.end(), list.getNumberOfNodes());
    m_numberOfNodes += list.getNumberOfNodes();
  }
}

void MeshIdNodeList::append(MeshIdNodeList list) {
  for (auto &it : list.m_data) {
    m_data.emplace_back(it.m_start, it.m_end, it.m_numberOfNodes);
    m_numberOfNodes += it.m_numberOfNodes;
  }
}

auto MeshIdNodeList::getDegreesOfFreedom() -> std::vector<DegreeOfFreedom *> {
  
  std::vector<DegreeOfFreedom *> dofs;
  dofs.reserve(m_numberOfNodes * 3);
  for (auto it = this->begin(); it != this->end();++it) {
    auto temp = it->getDegreesOfFreedom();
    dofs.insert(dofs.end(), temp.begin(), temp.end());
  }
  return dofs;
}

auto MeshIdNodeList::getNodeVector() -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nn;
  nn.reserve(m_numberOfNodes);
  for (auto it = this->begin(); it != this->end(); ++it) {
    nn.push_back(&(*it));
  }
  return nn;
}

auto MeshIdNodeList::operator[](indexType indx) -> GenericNodes & {
  if (indx >= m_numberOfNodes || indx < 0) {
    throw std::runtime_error("In MeshIdNodeList::operator[], accessed index out of bounds");
  } else {
    indexType pos = 0;
    indexType tindex = indx;
    while (tindex >= m_data[pos].m_numberOfNodes) {
      ++pos;
      tindex -= m_data[pos].m_numberOfNodes;
    }
    return *(m_data[pos].m_start + tindex);
  }

}

auto MeshIdNodeList::begin() -> iterator { return iterator(m_data.begin(),m_data.end()); }

auto MeshIdNodeList::end() -> iterator { return iterator(m_data.end(),m_data.end()); }

void MeshIdNodeList::add(std::vector<GenericNodes>::iterator startIter,
                         std::vector<GenericNodes>::iterator endIter,
                         indexType numberOfNodes) {

  m_data.emplace_back(startIter, endIter, numberOfNodes);
  m_numberOfNodes += numberOfNodes;
}

MeshIdNodeList::set::set(nodeiter start, nodeiter end, indexType numberOfNodes)
    : m_start(start), m_end(end), m_numberOfNodes(numberOfNodes) {}

auto MeshIdNodeList::set::begin() -> nodeiter { return m_start; }

auto MeshIdNodeList::set::end() -> nodeiter { return m_end; }

MeshIdNodeList::iterator::iterator(std::vector<set>::iterator start,
                                   std::vector<set>::iterator end)
    : m_outerIt(start), m_outerEnd(end) {
      update();
    }

MeshIdNodeList::iterator &MeshIdNodeList::iterator::operator++() {
  ++m_innerIt;
  if (m_innerIt == m_innerEnd) {
    ++m_outerIt;
    this->update();
  }

  return *this;
}

auto MeshIdNodeList::iterator::operator*() -> GenericNodes & {
  return *m_innerIt;
}

auto MeshIdNodeList::iterator::operator->() -> GenericNodes * {
  return &(*m_innerIt);
}


void MeshIdNodeList::iterator::update() {
  while (m_outerIt != m_outerEnd) {
    m_innerIt = (*m_outerIt).begin();
    m_innerEnd = (*m_outerIt).end();
    if (m_innerIt == m_innerEnd)
      ++m_outerIt;
    else
      break;
  }
}



bool operator==(const MeshIdNodeList::iterator &lhs,
                const MeshIdNodeList::iterator &rhs) {
  bool lhsEnd = lhs.m_outerIt == lhs.m_outerEnd;
  bool rhsEnd = rhs.m_outerIt == rhs.m_outerEnd;
  if (lhsEnd && rhsEnd)
    return true;
  if (lhsEnd != rhsEnd)
    return false;

  return (lhs.m_outerIt == rhs.m_outerIt) && (rhs.m_innerIt == rhs.m_innerIt);
  return false;
}

bool operator!=(const MeshIdNodeList::iterator &lhs,
                const MeshIdNodeList::iterator &rhs) {
  return !(lhs == rhs);
  ;
}

} /* namespace HierAMuS */
