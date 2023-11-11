// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once
#include <datatypes.h>



#include <vector>
#include "Nodetypes.h"
#include <ostream>

#include "DegreeOfFreedom.h"
#include <array>

namespace HierAMuS {

class GenericNodes {
public:
  GenericNodes(indexType startDofId, indexType nodeId);
  ~GenericNodes();
  void print(std::ostream &out);

  void setId(indexType nodeId) { this->m_nodeId = nodeId; };
  auto getId() -> indexType { return this->m_nodeId; };

  void initEqIds() {};

  void setDofStatus(indexType dof, dofStatus status);
  void setBoundaryCondition(indexType dof);
  void unsetBoundaryCondition(indexType dof);
  void setNodeType(NodeTypes type) { this->m_type = type; };

  auto getDegreesOfFreedom() -> std::vector<DegreeOfFreedom *>;
  void addDofsToVector(std::vector<DegreeOfFreedom *> &Dofs);

  auto getDegreeOfFreedom(indexType number) -> DegreeOfFreedom &;

  template<typename OStream>
  friend OStream & operator<<(OStream & os, const GenericNodes& c)
  { 
    std::string fmtString;
    fmt::format_to(std::ostream_iterator<char>(os), 
            "   Node id:  {:>12}, \n{}\n{}\n{}", 
            c.m_nodeId,c.m_dofs[0],c.m_dofs[1],c.m_dofs[2]);
   
    return os; 
  }
protected:
  NodeTypes m_type;
  indexType m_nodeId;
  std::array<DegreeOfFreedom,3> m_dofs;
  
};

} /* namespace HierAMuS */
