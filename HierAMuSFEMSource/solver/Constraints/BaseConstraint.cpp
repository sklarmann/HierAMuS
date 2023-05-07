// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "BaseConstraint.h"
#include "control/BinaryWrite.h"

namespace HierAMuS {

BaseConstraint::BaseConstraint() {}

BaseConstraint::BaseConstraint(const BaseConstraint &other) : m_id(other.m_id){}

void BaseConstraint::setId(indexType id) { m_id = id; }

auto BaseConstraint::getId() -> indexType { return m_id; }

void BaseConstraint::toFile(std::ofstream &out)
{ writeScalar(out, m_id); }

void BaseConstraint::fromFile(std::ifstream &in)
{ readScalar(in, m_id); }

} // namespace HierAMuS