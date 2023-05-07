// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "datatypes.h"

#include <equations/DofStatus.h>
#include <forwarddeclaration.h>

#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"

namespace HierAMuS {

class DegreeOfFreedom {
public:
  DegreeOfFreedom();
  ~DegreeOfFreedom() = default;

  void setId(indexType id);
  void setEqId(indexType equationId);
  void setStatus(dofStatus status);
  auto getId() const -> indexType;
  auto getEqId() const -> indexType;
  auto getStatus() const -> dofStatus;
  void print(PointerCollection &pointers);
  void setConstraintId(indexType id);
  auto getConstraintId() -> indexType;

  template<typename OStream>
  friend OStream & operator<<(OStream & os, const DegreeOfFreedom& c)
  { 
    std::string fmtString;
    if(c.m_status == dofStatus::active) {
      fmt::format_to(std::ostream_iterator<char>(os), 
            "Degree of Freedom:  {:>8},  Equation id:  {:>8},  Status:  active", 
            c.m_id,c.m_eqid);
    } else if(c.m_status == dofStatus::inactive){
      fmt::format_to(std::ostream_iterator<char>(os), 
            "Degree of Freedom:  {:>8},  Equation id:  {:>8},  Status:  inactive", 
            c.m_id,c.m_eqid);
    } else if(c.m_status == dofStatus::constraint) {
      fmt::format_to(std::ostream_iterator<char>(os), 
            "Degree of Freedom:  {:>8},  Equation id:  {:>8},  Status:  Dof is constraint by constraint  {:>8}", 
            c.m_id,c.m_eqid,c.m_constraintId);

    } else {
      fmt::format_to(std::ostream_iterator<char>(os), 
            "Degree of Freedom:  {:>8},  Equation id:  {:>8},  Status:  unknown", 
            c.m_id,c.m_eqid);
    }
    return os; 
  }

private:
  indexType m_id;
  indexType m_eqid;
  dofStatus m_status;
  indexType m_constraintId;
};


} // namespace HierAMuS
