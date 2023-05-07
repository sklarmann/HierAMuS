// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include <equations/NodeSet.h>

#include <equations/GenericNodes.h>
#include <equations/EquationHandler.h>

#include <ostream>
#include <pointercollection/pointercollection.h>
#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <control/FEMExceptions.h>

#include <iomanip>
#include <sstream>

namespace HierAMuS {



  
  NodeSet::NodeSet(){
    this->type     			= undef;
    this->numberOfNodes 	= 0;
    this->meshId 			= 0;
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
          nodeSetSetTypeException except;
          throw except;
        }
    }
  }

  void NodeSet::setMeshId(indexType meshId) { this->meshId = meshId; }

  indexType NodeSet::getMeshId() {
    // TODO: hier return-Anweisung eingeben
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
      nodeSetNumberOfNodesException except;
      throw except;
    }
  }

  
  void NodeSet::print(PointerCollection &pointers)
  {
    auto Logger = pointers.getSPDLogger();
    Logger.debug(*this);


    std::vector<GenericNodes*> nodes;
    pointers.getEquationHandler()->getNodes(nodes, (*this));
    for (auto & node : nodes) {
      Logger.debug(*node);
    }
  }


  
  void NodeSet::getNodes(
      PointerCollection &pointers,
      std::vector<GenericNodes*> &Nodes) {

    pointers.getEquationHandler()->getNodes(Nodes, *this);
  }

  auto NodeSet::getNodes(PointerCollection &pointers)
      -> std::vector<GenericNodes *> {
    std::vector<GenericNodes *> Nodes;
    pointers.getEquationHandler()->getNodes(Nodes, *this);
    return Nodes;
  }

  
  std::vector<DegreeOfFreedom*> HierAMuS::NodeSet::getDegreesOfFreedom(
      PointerCollection &pointers) {

    std::vector<DegreeOfFreedom*> ret,temp;
    std::vector<GenericNodes*> tnodes;
    this->getNodes(pointers,tnodes);
    for(auto & tnode : tnodes){
      temp.clear();
      tnode->getDegreesOfFreedom(pointers,temp);
      ret.insert(ret.end(), temp.begin(),temp.end());
    }


    return ret;
  }



} /* namespace HierAMuS */





