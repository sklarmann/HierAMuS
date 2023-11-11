// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include "NodeSet.h"

#include "GenericNodes.h"
#include "EquationHandler.h"

#include <ostream>


#include <iomanip>
#include <sstream>

namespace HierAMuS {



  
  NodeSet::NodeSet()
    : type(NodeTypes::undef), meshId(0), numberOfNodes(0), nodeStorageId(-1) {
  }

  NodeSet::NodeSet(NodeTypes typei, indexType meshIdi, indexType numberOfNodesi)
      : type(typei), meshId(meshIdi), numberOfNodes(numberOfNodesi),
        nodeStorageId(-1) {
  }


  
  NodeSet::~NodeSet(){

  }


  /**
   * @brief Set the node type of the set.
   * @param type Node type of the set
   */
  
  void NodeSet::setType(NodeTypes type){
    switch(this->type){
      case undef:
        this->type = type;
        break;
      default:
        if(this->type!=type){
          throw std::runtime_error("Switching node type not allowed");
        }
    }
  }

  void NodeSet::setMeshId(indexType meshId) { this->meshId = meshId; }

  indexType NodeSet::getMeshId() {
    return this->meshId; 
  }


  /**
   * @brief Set the number of nodes in the set.
   * @param numberOfNodes Number of nodes in the set.
   */
  
  void NodeSet::setNumberOfNodes(indexType numberOfNodes){

    if(this->numberOfNodes==0||this->numberOfNodes==numberOfNodes){
      this->numberOfNodes = numberOfNodes;
    }else{
      throw std::runtime_error("Changing the number of nodes in the set is not allowed!");
    }
  }


} /* namespace HierAMuS */





