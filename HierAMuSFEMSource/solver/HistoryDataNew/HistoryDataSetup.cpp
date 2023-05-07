// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "HistoryDataSetup.h"

namespace HierAMuS {
HistoryDataSetup::HistoryDataSetup()
    : m_elementId(-1), m_numberOfIntegrationPoints(0),
      m_numberOfHistoryVariablesUpdateElement(0),
      m_numberOfHistoryVariablesConstElement(0),
      m_numberOfHistoryVariablesUpdateMaterial(0),
      m_numberOfHistoryVariablesConstMaterial(0) {}

void HistoryDataSetup::setNumberOfIntegrationPoints(
    indexType numberOfIntegrationPoints) {
  m_numberOfIntegrationPoints = numberOfIntegrationPoints;
}

void HistoryDataSetup::setElementId(indexType id) { m_elementId = id; }

void HistoryDataSetup::setNHistPerIntPointElement(indexType update,
                                                  indexType constant) {
  if (m_numberOfHistoryVariablesUpdateElement < update)
    m_numberOfHistoryVariablesUpdateElement = update;
  if (m_numberOfHistoryVariablesConstElement < constant)
    m_numberOfHistoryVariablesConstElement = constant;
}

void HistoryDataSetup::setNHistPerIntPointMaterial(indexType update,
                                                   indexType constant) {
  if (m_numberOfHistoryVariablesUpdateMaterial < update)
    m_numberOfHistoryVariablesUpdateMaterial = update;
  if (m_numberOfHistoryVariablesConstMaterial < constant)
    m_numberOfHistoryVariablesConstMaterial = constant;
}

auto HistoryDataSetup::getTotalHistoryElementUpdate() -> indexType {
  return m_numberOfIntegrationPoints * m_numberOfHistoryVariablesUpdateElement;
}

auto HistoryDataSetup::getTotalHistoryElementConst() -> indexType {
  return m_numberOfIntegrationPoints * m_numberOfHistoryVariablesConstElement;
}

auto HistoryDataSetup::getTotalHistoryMaterialUpdate() -> indexType {
  return m_numberOfIntegrationPoints * m_numberOfHistoryVariablesUpdateMaterial;
}

auto HistoryDataSetup::getTotalHistoryMaterialConst() -> indexType {
  return m_numberOfIntegrationPoints * m_numberOfHistoryVariablesConstMaterial;
}

auto HistoryDataSetup::getElementId() -> indexType { return m_elementId; }

} // namespace HierAMuS