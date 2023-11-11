// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include "NodeSetList.h"
#include "NodeSet.h"

#include <ostream>

#include <iomanip>
#include <sstream>

namespace HierAMuS {

NodeSetList::NodeSetList(indexType numberOfSets,
                         std::vector<NodeSet>::iterator &startIterator,
                         std::vector<NodeSet>::iterator &endIterator) : numberOfNodeSets(numberOfSets), startIterator(startIterator), endIterator(endIterator) {}

NodeSetList::NodeSetList(indexType numberOfSets,
                         std::vector<NodeSet>::iterator &&startIterator,
                         std::vector<NodeSet>::iterator &&endIterator) : numberOfNodeSets(numberOfSets), startIterator(startIterator), endIterator(endIterator) {}


  
  NodeSetList::~NodeSetList(){

  }

auto NodeSetList::operator[](indexType indx) -> NodeSet & {
  // TODO: hier return-Anweisung eingeben
  if (indx < numberOfNodeSets) {
    return *(startIterator + indx);
  }
  throw std::runtime_error("In NodeSetList index is out of bounds!");
}

auto NodeSetList::hasSetMeshId(indexType meshId) -> bool { 
  
  for (auto it = this->startIterator; it != this->endIterator; ++it) {
    if (it->getMeshId() == meshId)
      return true;
  }
  return false;
}

auto NodeSetList::getSetMeshId(indexType meshId) -> NodeSet & {
  for (auto it = this->startIterator; it != this->endIterator; ++it) {
    if (it->getMeshId() == meshId)
      return *it;
  }

  throw std::runtime_error("In NodeSetList index is out of bounds!");
  //return NodeSet(NodeTypes::undef, meshId, 0);
}

auto NodeSetList::getNumberOfSets() -> indexType { return this->numberOfNodeSets; }

auto NodeSetList::begin() -> std::vector<NodeSet>::iterator {
  return startIterator;
}

auto NodeSetList::end() -> std::vector<NodeSet>::iterator {
  return endIterator;
}


} /* namespace HierAMuS */





