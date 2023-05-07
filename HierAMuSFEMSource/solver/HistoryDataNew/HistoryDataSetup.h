// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "datatypes.h"

namespace HierAMuS {

class HistoryDataSetup {

public:
  HistoryDataSetup();
  ~HistoryDataSetup() = default;

  void setNumberOfIntegrationPoints(indexType numberOfIntegrationPoints);
  void setElementId(indexType id);
  void setNHistPerIntPointElement(indexType update, indexType constant);
  void setNHistPerIntPointMaterial(indexType update, indexType constant);

  auto getTotalHistoryElementUpdate() -> indexType;
  auto getTotalHistoryElementConst() -> indexType;
    
  auto getTotalHistoryMaterialUpdate() -> indexType;
  auto getTotalHistoryMaterialConst() -> indexType;

  auto getElementId() -> indexType;


private:
  indexType m_elementId;

  indexType m_numberOfIntegrationPoints; // Number of integration points per element
  indexType
      m_numberOfHistoryVariablesUpdateElement; // Number of history variables
                                               // per integration point
  indexType
      m_numberOfHistoryVariablesConstElement; // Number of history variables per
                                              // integration point
  indexType
      m_numberOfHistoryVariablesUpdateMaterial; // Number of history variables
                                               // per integration point
  indexType
      m_numberOfHistoryVariablesConstMaterial; // Number of history variables per
                                              // integration point
};


} // namespace HierAMuS