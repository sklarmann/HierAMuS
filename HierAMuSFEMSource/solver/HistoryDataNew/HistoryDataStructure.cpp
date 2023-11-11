// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "HistoryDataStructure.h"


namespace HierAMuS {



HistoryDataStructure::HistoryDataStructure(
  std::vector<std::pair<indexType, indexType>> constHistory,
  std::vector<std::pair<indexType, indexType>> updateHistory)
    : m_constHistory(std::move(constHistory)), m_updateHistory(std::move(updateHistory))
{
  update_number_of_const_values();
  update_number_of_update_values();
}

auto HistoryDataStructure::getNumberOfUpdateValues() const -> indexType {
  return m_numberOfUpdateValues;
}

auto HistoryDataStructure::getNumberOfConstValues() const -> indexType {
  return m_numberOfConstValues;
}

auto HistoryDataStructure::getConstStructure() const
    -> const std::vector<std::pair<indexType, indexType>>& {
  return m_constHistory;
}

auto HistoryDataStructure::getUpdateStructure() const
    -> const std::vector<std::pair<indexType, indexType>> & {
  return m_updateHistory;
}

void HistoryDataStructure::update_number_of_update_values() {
  m_numberOfUpdateValues = 0;
  for (auto &i : m_updateHistory) {
    m_numberOfUpdateValues += i.first * i.second;
  }
}

void HistoryDataStructure::update_number_of_const_values() {
  m_numberOfConstValues = 0;
  for (auto &i : m_constHistory) {
    m_numberOfConstValues += i.first * i.second;
  }
}

} // namespace HierAMuS
