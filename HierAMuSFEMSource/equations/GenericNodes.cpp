// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include <equations/GenericNodes.h>

#include <equations/DegreeOfFreedom.h>
#include <equations/DofStatus.h>
#include <equations/EquationHandler.h>
#include <equations/Nodetypes.h>
#include <iostream>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>

#include <pointercollection/pointercollection.h>

#include <iomanip>

namespace HierAMuS {

GenericNodes::GenericNodes(indexType &startDofId) {
  this->m_type = NodeTypes::undef;
  for (auto &i : m_dofs) {
    i.setId(startDofId);
    ++startDofId;
  }
}

GenericNodes::~GenericNodes() = default;

void GenericNodes::print(std::ostream &out) {
  out << "Generic Node print! Id: " << this->m_nodeId << std::endl;
}

void GenericNodes::setDofStatus(indexType dof, dofStatus status) {
  m_dofs[dof].setStatus(status);
};

void GenericNodes::setBoundaryCondition(PointerCollection &pointers,
                                        indexType dof) {
  if (dof > 2) {
    // TODO throw exception
    std::cout << "dof exceeds number of degrees of freedom" << std::endl;
  }
  m_dofs[dof].setStatus(dofStatus::inactive);
}

void GenericNodes::unsetBoundaryCondition(PointerCollection &pointers,
                                          indexType dof) {
  m_dofs[dof].setStatus(dofStatus::active);
}

void GenericNodes::getDegreesOfFreedom(PointerCollection &pointers,
                                       std::vector<DegreeOfFreedom *> &Dofs) {
  Dofs = {&m_dofs[0], &m_dofs[1], &m_dofs[2]};
}

auto GenericNodes::getDegreesOfFreedom(PointerCollection &pointers)
    -> std::vector<DegreeOfFreedom *> {

  std::vector<DegreeOfFreedom *> Dofs({&m_dofs[0], &m_dofs[1], &m_dofs[2]});
  return Dofs;
}

auto GenericNodes::getDegreeOfFreedom(indexType number) -> DegreeOfFreedom & {
  if (number > 2) {
    throw std::runtime_error(
        "In Node: dof number exceeds number of degrees of freedom");
  }
  return this->m_dofs[number];
}

void GenericNodes::print(PointerCollection &pointers) {
  pointers.getSPDLogger().debug(*this);
  
}

} /* namespace HierAMuS */
