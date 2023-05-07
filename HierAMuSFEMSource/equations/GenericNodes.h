// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once
#include <datatypes.h>

#include <forwarddeclaration.h>

#include <vector>
#include <equations/Nodetypes.h>
#include <ostream>

#include "equations/DegreeOfFreedom.h"
#include <array>

namespace HierAMuS {

class GenericNodes {
public:
  GenericNodes(indexType &startDofId);
  ~GenericNodes();
  void print(std::ostream &out);

  void setId(indexType nodeId) { this->m_nodeId = nodeId; };
  auto getId() -> indexType { return this->m_nodeId; };

  void initEqIds() {};

  void setDofStatus(indexType dof, dofStatus status);
  void setBoundaryCondition(PointerCollection &pointers, indexType dof);
  void unsetBoundaryCondition(PointerCollection &pointers, indexType dof);
  void setNodeType(NodeTypes type) { this->m_type = type; };
  void getDegreesOfFreedom(PointerCollection &pointers,
                           std::vector<DegreeOfFreedom *> &Dofs);
  auto getDegreesOfFreedom(PointerCollection &pointers) -> std::vector<DegreeOfFreedom *>;

  auto getDegreeOfFreedom(indexType number) -> DegreeOfFreedom &;


  void print(PointerCollection &pointers);


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
  indexType m_dofStorageId;
  std::array<DegreeOfFreedom,3> m_dofs;
  
};

} /* namespace HierAMuS */
