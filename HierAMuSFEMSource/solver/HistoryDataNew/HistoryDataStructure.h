// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "datatypes.h"

namespace HierAMuS {

class HistoryDataStructure {

public:
  HistoryDataStructure(
    std::vector<std::pair<indexType, indexType>> constHistory,
    std::vector<std::pair<indexType, indexType>> updateHistory);
  ~HistoryDataStructure() = default;

  auto getNumberOfUpdateValues() const -> indexType;
  auto getNumberOfConstValues() const -> indexType;

  auto getConstStructure() const
      -> const std::vector<std::pair<indexType, indexType>>&;
  auto getUpdateStructure() const
      -> const std::vector<std::pair<indexType, indexType>> &;
  
private:
  std::vector<std::pair<indexType, indexType>> m_constHistory, m_updateHistory;
  indexType m_numberOfConstValues;
  indexType m_numberOfUpdateValues;
};


} // namespace HierAMuS