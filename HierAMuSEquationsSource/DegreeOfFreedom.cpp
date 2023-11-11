// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "DegreeOfFreedom.h"
#include <iomanip>


namespace HierAMuS {


  
  DegreeOfFreedom::DegreeOfFreedom()
    : m_id(-1), m_eqid(-1), m_status(dofStatus::active), m_constraintId(-1) {}



  void DegreeOfFreedom::setConstraintId(indexType id) { m_constraintId = id; }
  auto DegreeOfFreedom::getConstraintId() -> indexType { return m_constraintId; }
  void DegreeOfFreedom::setId(indexType id) { m_id = id; }
  void DegreeOfFreedom::setEqId(indexType equationId) {
    m_eqid = equationId; }
  void DegreeOfFreedom::setStatus(dofStatus status) { m_status = status; }
  }

