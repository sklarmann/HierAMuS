// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "datatypes.h"
#include "types/MatrixTypes.h"
#include "HistoryDataSetup.h"
#include "HistoryDataIterator.h"

namespace HierAMuS {

class HistoryDataManager {

public:
  HistoryDataManager() = default;
  ~HistoryDataManager() = default;

  void setNumberOfElements(indexType numberOfElements);
  void addHistoryDataElement(HistoryDataSetup &historyDataStructure);

  void initHistoryData();

  auto getHistoryDataIterator(indexType elementId) -> HistoryDataIterator;

  void update();

  void toFile(std::ofstream &out);
  void fromFile(std::ifstream &in);

private:
  // Update element history data
  std::vector<indexType> m_historyElementUpdatePositions;
  Types::VectorX<prec> m_historyElementUpdate_old;
  Types::VectorX<prec> m_historyElementUpdate_new;

  std::vector<indexType> m_historyElementConstPositions;
  Types::VectorX<prec> m_historyElementConst;

  // Update material history data
  std::vector<indexType> m_historyMaterialUpdatePositions;
  Types::VectorX<prec> m_historyMaterialUpdate_old;
  Types::VectorX<prec> m_historyMaterialUpdate_new;

  std::vector<indexType> m_historyMaterialConstPositions;
  Types::VectorX<prec> m_historyMaterialConst;



};


} // namespace HierAMuS