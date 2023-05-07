// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#ifndef NODESETSTORAGE_H_
#define NODESETSTORAGE_H_

#include <NodeSet.h>
#include <vector>
#include <DofHandler.h>

namespace HierAMuS {

  
  class NodeSetStorage {
    public:
      NodeSetStorage(DofHandler *dofHandler);
      ~NodeSetStorage();
      void requestNodeSets(indexType &first, const unsigned char &num);
      NodeSet &getNodeSet(const indexType &num);
      void setDofHandler(DofHandler *dofHandlerin) {this->dofHandler = dofHandlerin;};
    private:
      std::vector<NodeSet> nodeSets;
      DofHandler *dofHandler;
  };


  
  NodeSetStorage::NodeSetStorage(DofHandler *dofHandler){
    this -> nodeSets.reserve(10000);
    this -> dofHandler = dofHandler;
  }


  
  NodeSetStorage::~NodeSetStorage(){
    nodeSets.clear();
  }

  
  void NodeSetStorage::requestNodeSets(indexType &first, const unsigned char &num){
    first = this->nodeSets.size();

    for(auto i=0;i<num;++i){
      this->nodeSets.emplace_back(&this->dofHandler);
    }
  }

  
  NodeSet &NodeSetStorage::getNodeSet(const indexType &num){
    if(num >= this->nodeSets.size()) {
      return 0;
    }
    return this->nodeSets.at(num);
  }



} /* namespace HierAMuS */



#endif /* NODESETSTORAGE_H_ */
