// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include <geometry/Base.h>

#include <control/HandlingStructs.h>
#include <pointercollection/pointercollection.h>

#include "geometry/GeometryData.h"

#include <equations/EquationHandler.h>
#include <equations/GenericNodes.h>
#include <equations/NodeSet.h>

#include <iomanip>

namespace HierAMuS::Geometry {

Base::Base() {
  this->numberOfNodeSets = 0;
  this->nodeSetId = 0;
  this->nodeSetInit = false;
}

Base::~Base() = default;

auto Base::getType() -> const GeometryTypes & { return HierAMuS::Geometry::Base::type; };
auto Base::getGroupType() -> const GeometryTypes & { return HierAMuS::Geometry::Base::type; };

void Base::setNodeSet(PointerCollection &pointers, indexType meshID,
                      indexType numberOfNodes, NodeTypes type) {
  this->nodeSetId = pointers.getEquationHandler()->requestNodeSetSetup(
      this->nodeSetId, meshID, this->nodeSetInit);
  this->numberOfNodeSets =
      pointers.getEquationHandler()->getNumberOfNodeSets(this->nodeSetId);

  NodeSet *tempSet =
      pointers.getEquationHandler()->getSetMeshId(this->nodeSetId, meshID);
  tempSet->setNumberOfNodes(numberOfNodes);
  tempSet->setType(type);
}

void Base::setVerts(GeometryData &geoData, std::vector<indexType> &vertsIn){}

// void Base::print(PointerCollection &pointers) {
//   OutputHandler &Log = pointers.getInfoData()->Log;
//   Log.debug()
//       << "Generic geometry element call" << std::left << std::setw(20)
//       << "\tId:" << std::right << std::setw(10) << this->id << "\n"
//       << std::left << std::setw(20) << "\tNumber of nodesets:" << std::right
//       << std::setw(10) << this->numberOfNodeSets << std::endl;
// }

void Base::printEqInfo(PointerCollection &pointers) {

  auto &Logger = pointers.getSPDLogger();

  Logger.debug("   Id: {:>12}    Number of nodesets:  {:>12}",
               this->id, this->numberOfNodeSets);

  if (this->numberOfNodeSets > 0) {
    std::vector<NodeSet *> sets;
    pointers.getEquationHandler()->getSets(sets, this->nodeSetId,
                                           this->numberOfNodeSets);
    for (auto & set : sets) {
      set->print(pointers);
    }
  }
}

auto Base::getSet(ptrCol &pointers, indexType setNumber) -> NodeSet * {
  if (setNumber >= this->numberOfNodeSets) {
    // TODO throw exceptions
  }
  return pointers.getEquationHandler()->getSet(this->nodeSetId, setNumber);
}

void Base::getSets(ptrCol &pointers, std::vector<NodeSet *> &sets) {
  sets.clear();
  for (auto i = 0; i < this->numberOfNodeSets; ++i) {
    sets.push_back(this->getSet(pointers, i));
  }
}

void Base::setAllNodeBoundaryConditionMeshId(ptrCol &pointers, indexType meshId,
                                             indexType dof) {
  if (this->nodeSetInit) {
    auto temp = pointers.getEquationHandler();
    NodeSet *tempNodeSet;
    tempNodeSet = temp->getSetMeshId(this->nodeSetId, meshId);

    if (tempNodeSet != nullptr) {
      GenericNodes *tempNode;
      for (auto i = 0; i < tempNodeSet->getNumberOfNodes(); ++i) {
        tempNode = temp->getNode(*tempNodeSet, i);
        tempNode->setBoundaryCondition(pointers, dof);
      }
    }
  }
}

void Base::getNodesOfSet(ptrCol &pointers,
                         std::vector<GenericNodes *> &nodeVector,
                         indexType meshId) {

  nodeVector.clear();
  if (this->nodeSetInit) {
    NodeSet *tempSet;
    tempSet = pointers.getEquationHandler()->getSetMeshId(
        this->nodeSetId, meshId); // ->getNodes(nodeVector);
    if (tempSet != nullptr) {
      pointers.getEquationHandler()->getNodes(nodeVector, *tempSet);
}
  }
}

auto Base::getNodesOfSet(ptrCol &pointers,
                                                indexType meshId) -> std::vector<GenericNodes *> {
  std::vector<GenericNodes *> nodeVector;
  if (this->nodeSetInit) {
    NodeSet *tempSet;
    tempSet = pointers.getEquationHandler()->getSetMeshId(
        this->nodeSetId, meshId); // ->getNodes(nodeVector);
    if (tempSet != nullptr) {
      nodeVector = tempSet->getNodes(pointers);
    }
  }
  return nodeVector;
}


void Base::getAllEquationsIds(ptrCol &pointers,
                              std::vector<DegreeOfFreedom *> &Dofs) {
  Dofs.clear();

  auto tempEq = pointers.getEquationHandler();
  for (auto i = 0; i < this->numberOfNodeSets; ++i) {
    unsigned char numNodes =
        tempEq->getSet(this->nodeSetId, i)->getNumberOfNodes();
    for (unsigned char j = 0; j < numNodes; ++j) {
      std::vector<DegreeOfFreedom *> tempDofs;
      tempEq->getNode((*(tempEq->getSet(this->nodeSetId, i))), j)
          ->getDegreesOfFreedom(pointers, tempDofs);
      for (auto & tempDof : tempDofs) {
        Dofs.push_back(tempDof);
      }
    }
  }
}

void Base::getNodeEquationIds(ptrCol &pointers,
                              std::vector<DegreeOfFreedom *> &Dofs,
                              indexType meshId, indexType nodeNumber) {
  NodeSet *tempSet = this->getSetMeshId(pointers, meshId);
  if (tempSet == nullptr) {
    return;
}
  auto numnodes = tempSet->getNumberOfNodes();
  if (nodeNumber <= numnodes) {
    GenericNodes *tempnode =
        pointers.getEquationHandler()->getNode(*tempSet, nodeNumber);
    tempnode->getDegreesOfFreedom(pointers, Dofs);
  }
}

auto Base::getSetMeshId(ptrCol &pointers, indexType meshId) -> NodeSet * {
  NodeSet *ret = nullptr;
  if (this->numberOfNodeSets > 0) {
    ret = pointers.getEquationHandler()->getSetMeshId(this->nodeSetId, meshId);
  }
  return ret;
}

auto Base::checkSetWithMeshId(ptrCol &pointers, indexType meshId) -> bool {
  if (this->numberOfNodeSets > 0) {
    bool ret = false;
    unsigned char i = 0;
    NodeSet *temp;
    while (i < this->numberOfNodeSets && !ret) {
      temp = pointers.getEquationHandler()->getSet(this->nodeSetId, i);
      if (meshId == temp->getMeshId()) {
        ret = true;
}
      ++i;
    }

    return ret;
  }     return false;

}

void Base::getNodes(Base::ptrCol &pointers,
                    std::vector<GenericNodes *> &nodeVector, indexType meshId) {
}

void Base::getNodesInternal(Base::ptrCol &pointers,
                            std::vector<GenericNodes *> &nodeVector,
                            indexType meshId) {}

auto Base::getIntegrationPoints(ptrCol &pointers, indexType elementId) -> IntegrationPoints {
  auto temp = HierAMuS::Geometry::Base::ptrCol::getIntegrationPoints(-1);
  return temp;
}

const HierAMuS::Geometry::GeometryTypes Base::type =
    HierAMuS::Geometry::GeometryTypes::Generic;

} // namespace HierAMuS
