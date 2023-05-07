// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "HistoryDataStructureCombinded.h"
#include "HistoryDataStructure.h"

namespace HierAMuS {
HistoryDataStructureCombinded::HistoryDataStructureCombinded(
    const HistoryDataStructure &elementHistory, const HistoryDataStructure &materialHistory)
    : m_elementHistoryStructure(elementHistory),
      m_materialHistoryStructure(materialHistory) {}
auto HistoryDataStructureCombinded::getElementHistoryStructure() const
    -> const HistoryDataStructure & {
  return m_elementHistoryStructure;
}
auto HistoryDataStructureCombinded::getMaterialHistoryStructure() const
    -> const HistoryDataStructure & {
  return m_materialHistoryStructure;
}
} // namespace HierAMuS