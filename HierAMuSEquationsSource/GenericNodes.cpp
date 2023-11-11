// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "GenericNodes.h"

#include "DegreeOfFreedom.h"
#include "DofStatus.h"
#include "EquationHandler.h"
#include "Nodetypes.h"
#include <iostream>



#include <iomanip>

namespace HierAMuS {

GenericNodes::GenericNodes(indexType startDofId, indexType nodeId)
    : m_type(NodeTypes::undef), m_nodeId(nodeId),
      m_dofs({DegreeOfFreedom(startDofId), DegreeOfFreedom(startDofId + 1),
              DegreeOfFreedom(startDofId + 2)}) {}

GenericNodes::~GenericNodes() = default;

void GenericNodes::print(std::ostream &out) {
  out << "Generic Node print! Id: " << this->m_nodeId << std::endl;
}

void GenericNodes::setDofStatus(indexType dof, dofStatus status) {
  m_dofs[dof].setStatus(status);
};

void GenericNodes::setBoundaryCondition(indexType dof) {
  if (dof > 2) {
    // TODO throw exception
    std::cout << "dof exceeds number of degrees of freedom" << std::endl;
  }
  m_dofs[dof].setStatus(dofStatus::inactive);
}

void GenericNodes::unsetBoundaryCondition(indexType dof) {
  m_dofs[dof].setStatus(dofStatus::active);
}

auto GenericNodes::getDegreesOfFreedom() -> std::vector<DegreeOfFreedom *> {
  return {&m_dofs[0], &m_dofs[1], &m_dofs[2]};
}

void GenericNodes::addDofsToVector(std::vector<DegreeOfFreedom *> &Dofs) {
  std::vector<DegreeOfFreedom *> tt = {&m_dofs[0], &m_dofs[1], &m_dofs[2]};
  Dofs.insert(Dofs.end(), tt.begin(), tt.end());
}

auto GenericNodes::getDegreeOfFreedom(indexType number) -> DegreeOfFreedom & {
  if (number > 2) {
    throw std::runtime_error(
        "In Node: dof number exceeds number of degrees of freedom");
  }
  return this->m_dofs[number];
}

} /* namespace HierAMuS */
