// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <equations/DegreeOfFreedom.h>
#include <control/OutputHandler.h>
#include <iomanip>
#include <pointercollection/pointercollection.h>
#include <control/HandlingStructs.h>


namespace HierAMuS {


  
  DegreeOfFreedom::DegreeOfFreedom()
    : m_status(dofStatus::active), m_eqid(-1), m_id(-1), m_constraintId(-1) {}

  void DegreeOfFreedom::print(PointerCollection& pointers)
  {
    pointers.getSPDLogger().debug(*this);
  }

  void DegreeOfFreedom::setConstraintId(indexType id) { m_constraintId = id; }
  auto DegreeOfFreedom::getConstraintId() -> indexType { return m_constraintId; }
  auto DegreeOfFreedom::getId() const -> indexType { return m_id; }
  auto DegreeOfFreedom::getEqId() const -> indexType { return m_eqid; }
  auto DegreeOfFreedom::getStatus() const -> dofStatus { return this->m_status; }
  void DegreeOfFreedom::setId(indexType id) { m_id = id; }
  void DegreeOfFreedom::setEqId(indexType equationId) {
    m_eqid = equationId; }
  void DegreeOfFreedom::setStatus(dofStatus status) { m_status = status; }
  }

