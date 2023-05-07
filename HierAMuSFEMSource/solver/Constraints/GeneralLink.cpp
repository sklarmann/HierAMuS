// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "GeneralLink.h"

#include "equations/DegreeOfFreedom.h"
#include "equations/EquationHandler.h"

#include "control/BinaryWrite.h"

namespace HierAMuS {

GeneralLink::GeneralLink() {}

GeneralLink::GeneralLink(const GeneralLink &other)
    : BaseConstraint(other), m_masterDof(other.m_masterDof),
      m_slaveDof(other.m_slaveDof), m_factor(other.m_factor),
      m_difference(other.m_difference) {}

auto GeneralLink::getCopy() -> std::shared_ptr<BaseConstraint> {
  return std::make_shared<GeneralLink>(GeneralLink(*this));
}

auto GeneralLink::getDofs(PointerCollection &pointers)
    -> std::vector<DegreeOfFreedom *> {
  std::vector<DegreeOfFreedom *> Dofs(
      {&pointers.getEquationHandler()->getDegreeOfFreedom(m_masterDof),
       &pointers.getEquationHandler()->getDegreeOfFreedom(m_slaveDof)});

  return Dofs;
}

auto GeneralLink::getDofRelation(PointerCollection &pointers)
    -> DofRelation {

  DofRelation dofsr;
  auto eqHandler = pointers.getEquationHandler();
  dofsr.m_slaveDofs.push_back(&eqHandler->getDegreeOfFreedom(m_slaveDof));
  dofsr.m_masterDofs.push_back(&eqHandler->getDegreeOfFreedom(m_masterDof));
  
  return dofsr;
}

void GeneralLink::set(PointerCollection &pointers, indexType masterDof,
                      indexType slaveDof, prec factor, prec difference) {
  m_masterDof = masterDof;
  m_slaveDof = slaveDof;
  m_factor = factor;
  m_difference = difference;
  auto &Dof = pointers.getEquationHandler()->getDegreeOfFreedom(m_slaveDof);
  Dof.setStatus(dofStatus::constraint);
  Dof.setConstraintId(this->getId());
  auto &Dofm = pointers.getEquationHandler()->getDegreeOfFreedom(m_masterDof);
  Dofm.setConstraintId(this->getId());
  
 
}

auto GeneralLink::getSlaveDisplacement(PointerCollection &pointers)
    -> Types::SparseVector<prec, indexType> {

  auto eqHandler = pointers.getEquationHandler();
  Types::SparseVector<prec, indexType> slaveDisp(eqHandler->getNumberOfTotalEquations());

  auto mDof = eqHandler->getDegreeOfFreedom(m_masterDof);
  auto sDof = eqHandler->getDegreeOfFreedom(m_slaveDof);

  prec sSol = pointers.getSolutionState()->getSolution(sDof.getId());
  std::vector<DegreeOfFreedom *> Dofs;
  Dofs.push_back(&mDof);
  Types::VectorX<prec> sol =
      pointers.getSolutionState()->getNewtonSolution(Dofs);
  sSol += m_factor * sol[0];
  sSol += m_difference;
  slaveDisp.insert(sDof.getId()) = sSol;
  return slaveDisp;
}

auto GeneralLink::getSlaveNewtonSolution(PointerCollection &pointers)
    -> Types::SparseVector<prec, indexType> {
  auto eqHandler = pointers.getEquationHandler();
  Types::SparseVector<prec, indexType> slaveDisp(eqHandler->getNumberOfTotalEquations());

  auto mDof = eqHandler->getDegreeOfFreedom(m_masterDof);
  auto sDof = eqHandler->getDegreeOfFreedom(m_slaveDof);

  auto msol = pointers.getSolutionState()->getSolution(mDof.getId());
  auto ssol = pointers.getSolutionState()->getSolution(sDof.getId());
  prec dB = m_difference - ssol + m_factor * msol;

  slaveDisp.insert(sDof.getId()) = dB;

  return slaveDisp;
}

void GeneralLink::modifyEquationSystem(
    PointerCollection &pointers, Eigen::SparseMatrix<prec, 0, indexType> &Kaa,
    Eigen::SparseMatrix<prec, 0, indexType> &Kab,
    Eigen::SparseMatrix<prec, 0, indexType> &Kba,
    Eigen::SparseMatrix<prec, 0, indexType> &Kbb, Types::VectorX<prec> &Fa,
    Types::VectorX<prec> &Fb) {
  auto &DofM = pointers.getEquationHandler()->getDegreeOfFreedom(m_masterDof);

  if (DofM.getStatus() == dofStatus::active) {
    auto &DofS = pointers.getEquationHandler()->getDegreeOfFreedom(m_slaveDof);

    auto col = Types::SparseVector<prec, indexType>(Kab.col(DofS.getEqId()));
    auto row = Types::SparseVector<prec, indexType>(Kba.row(DofS.getEqId()));
    prec KBBval = Kbb.coeff(DofS.getEqId(), DofS.getEqId());

    auto msol = pointers.getSolutionState()->getSolution(DofM.getId());
    auto ssol = pointers.getSolutionState()->getSolution(DofS.getId());
    prec dB = m_difference - ssol + m_factor * msol;
    dB *= prec(-1);

    Fa += col * dB;
    KBBval *= m_factor;
    Fa(DofM.getEqId()) += KBBval * dB;
    KBBval *= m_factor;
    col *= m_factor;
    row *= m_factor;

    //col.coeffRef(DofM.getEqId()) += KBBval;
    pointers.getSolutionState()->insertRowInEqSystem(Kaa, row, DofM.getEqId());
    pointers.getSolutionState()->insertColumnInEqSystem(Kaa, col,
                                                        DofM.getEqId());
    Kaa.coeffRef(DofM.getEqId(), DofM.getEqId()) += KBBval;
    auto temp = (m_factor * Fb(DofS.getEqId()));
    Fa(DofM.getEqId()) += temp;
  }
}

auto GeneralLink::getAMatrix(PointerCollection &pointers)
    -> Types::SparseMatrix<prec, indexType>  {
  Types::SparseMatrix<prec, indexType> AMat;
  auto eqHandler = pointers.getEquationHandler();
  AMat.resize(eqHandler->getNumberOfInActiveEquations(),
              eqHandler->getNumberOfActiveEquations());
  auto Dofs = this->getDofRelation(pointers);
  AMat.coeffRef(Dofs.m_slaveDofs[0]->getEqId(),
                Dofs.m_masterDofs[0]->getEqId()) = m_factor;
  AMat.makeCompressed();
  return AMat;
}

void GeneralLink::setB(prec B)
{ m_difference = B; }

auto GeneralLink::getDB(PointerCollection &pointers)
    -> Types::SparseVector<prec, indexType> {

  auto eqHandler = pointers.getEquationHandler();
  Types::SparseVector<prec, indexType> dB;
  dB.resize(eqHandler->getNumberOfInActiveEquations());


  auto &DofM = pointers.getEquationHandler()->getDegreeOfFreedom(m_masterDof);
  
  auto &DofS = pointers.getEquationHandler()->getDegreeOfFreedom(m_slaveDof);
  auto msol = pointers.getSolutionState()->getSolution(DofM.getId());
  auto ssol = pointers.getSolutionState()->getSolution(DofS.getId());
  prec dBB = m_difference - ssol + m_factor * msol;
  //dBB *= prec(-1);

  dB.coeffRef(DofS.getEqId()) = dBB;
  

  return dB;
}

void GeneralLink::toFile(std::ofstream &out) {
  BaseConstraint::toFile(out);
  writeScalar(out, m_masterDof);
  writeScalar(out, m_slaveDof);
  writeScalar(out, m_factor);
  writeScalar(out, m_difference);
}

void GeneralLink::fromFile(std::ifstream &in) {
  BaseConstraint::fromFile(in);
  readScalar(in, m_masterDof);
  readScalar(in, m_slaveDof);
  readScalar(in, m_factor);
  readScalar(in, m_difference);
}

} // namespace HierAMuS