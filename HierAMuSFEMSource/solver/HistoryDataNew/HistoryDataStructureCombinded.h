// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace HierAMuS {

class HistoryDataStructure;
class HistoryDataStructureCombinded {

public:
  HistoryDataStructureCombinded(const HistoryDataStructure &elementHistory,
                                const HistoryDataStructure &materialHistory);


  auto getElementHistoryStructure() const -> const HistoryDataStructure &;
  auto getMaterialHistoryStructure() const -> const HistoryDataStructure &;

private:
  
  const HistoryDataStructure &m_elementHistoryStructure;
  const HistoryDataStructure &m_materialHistoryStructure;
};


} // namespace HierAMuS