// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "HistoryDataManager.h"
#include "control/BinaryWrite.h"


namespace HierAMuS {

void HistoryDataManager::setNumberOfElements(indexType numberOfElements) {
  // Element Data
  m_historyElementUpdatePositions.resize(numberOfElements + 1);
  for (auto &i : m_historyElementUpdatePositions) {
    i = 0;
  }
  m_historyElementConstPositions.resize(numberOfElements + 1);
  for (auto &i : m_historyElementConstPositions) {
    i = 0;
  }

  // Material Data
  m_historyMaterialConstPositions.resize(numberOfElements + 1);
  for (auto &i : m_historyMaterialConstPositions) {
    i = 0;
  }
  m_historyMaterialUpdatePositions.resize(numberOfElements + 1);
  for (auto &i : m_historyMaterialUpdatePositions) {
    i = 0;
  }
}

void HierAMuS::HistoryDataManager::addHistoryDataElement(
    HistoryDataSetup &historyDataStructure) {
  indexType elementId = historyDataStructure.getElementId();
  m_historyElementUpdatePositions[elementId + 1] =
      m_historyElementUpdatePositions[elementId] +
      historyDataStructure.getTotalHistoryElementUpdate();
  m_historyElementConstPositions[elementId + 1] =
      m_historyElementConstPositions[elementId] +
      historyDataStructure.getTotalHistoryElementConst();

  m_historyMaterialConstPositions[elementId + 1] =
      m_historyMaterialConstPositions[elementId] +
      historyDataStructure.getTotalHistoryMaterialConst();
  m_historyMaterialUpdatePositions[elementId + 1] =
      m_historyMaterialUpdatePositions[elementId] +
      historyDataStructure.getTotalHistoryMaterialUpdate();
}

void HistoryDataManager::initHistoryData() {
  indexType tothis = m_historyElementUpdatePositions.back();
  m_historyElementUpdate_new.resize(tothis);
  m_historyElementUpdate_old.resize(tothis);
  m_historyElementUpdate_new.setZero();
  m_historyElementUpdate_old.setZero();

  tothis = m_historyElementConstPositions.back();
  m_historyElementConst.resize(tothis);
  m_historyElementConst.setZero();

  tothis = m_historyMaterialUpdatePositions.back();
  m_historyMaterialUpdate_new.resize(tothis);
  m_historyMaterialUpdate_old.resize(tothis);
  m_historyMaterialUpdate_new.setZero();
  m_historyMaterialUpdate_old.setZero();

  tothis = m_historyMaterialConstPositions.back();
  m_historyMaterialConst.resize(tothis);
  m_historyMaterialConst.setZero();
}

auto HistoryDataManager::getHistoryDataIterator(indexType elementId)
    -> HistoryDataIterator {

  // Element Update History
  indexType elementUpdateStart = m_historyElementUpdatePositions[elementId];
  indexType elementUpdateEnd = m_historyElementUpdatePositions[elementId + 1];
  indexType elementUpdateLength = elementUpdateEnd - elementUpdateStart;
  Eigen::Map<Types::VectorX<prec>> elementUpdate_new(
      m_historyElementUpdate_new.data() + elementUpdateStart,
      elementUpdateLength);
  Eigen::Map<Types::VectorX<prec>> elementUpdate_old(
      m_historyElementUpdate_old.data() + elementUpdateStart,
      elementUpdateLength);

  // element const
  elementUpdateStart = m_historyElementConstPositions[elementId];
  elementUpdateEnd = m_historyElementConstPositions[elementId + 1];
  elementUpdateLength = elementUpdateEnd - elementUpdateStart;
  Eigen::Map<Types::VectorX<prec>> elementConst(
      m_historyElementConst.data() + elementUpdateStart, elementUpdateLength);

  // Material Update History
  elementUpdateStart = m_historyMaterialUpdatePositions[elementId];
  elementUpdateEnd = m_historyMaterialUpdatePositions[elementId + 1];
  elementUpdateLength = elementUpdateEnd - elementUpdateStart;
  Eigen::Map<Types::VectorX<prec>> materialUpdate_new(
      m_historyMaterialUpdate_new.data() + elementUpdateStart,
      elementUpdateLength);
  Eigen::Map<Types::VectorX<prec>> materialUpdate_old(
      m_historyMaterialUpdate_old.data() + elementUpdateStart,
      elementUpdateLength);

  // element const
  elementUpdateStart = m_historyMaterialConstPositions[elementId];
  elementUpdateEnd = m_historyMaterialConstPositions[elementId + 1];
  elementUpdateLength = elementUpdateEnd - elementUpdateStart;
  Eigen::Map<Types::VectorX<prec>> materialConst(
      m_historyMaterialConst.data() + elementUpdateStart, elementUpdateLength);

  HistoryDataIterator it(elementUpdate_old, elementUpdate_new, elementConst,
                         materialUpdate_new, materialUpdate_old, materialConst);

  return it;
}

void HistoryDataManager::update()
{
  m_historyElementUpdate_old = m_historyElementUpdate_new;
  m_historyMaterialUpdate_old = m_historyMaterialUpdate_new;
}

void HistoryDataManager::toFile(std::ofstream &out) {
  writeStdVector(out, m_historyElementUpdatePositions);
  writeEigenMatrix(out, m_historyElementUpdate_old);
  writeEigenMatrix(out, m_historyElementUpdate_new);
  writeStdVector(out, m_historyElementConstPositions);
  writeEigenMatrix(out, m_historyElementConst);
  writeStdVector(out, m_historyMaterialUpdatePositions);
  writeEigenMatrix(out, m_historyMaterialUpdate_old);
  writeEigenMatrix(out, m_historyMaterialUpdate_new);
  writeStdVector(out, m_historyMaterialConstPositions);
  writeEigenMatrix(out, m_historyMaterialConst);
}

void HistoryDataManager::fromFile(std::ifstream &in)
{
  readStdVector(in, m_historyElementUpdatePositions);
  readEigenMatrix(in, m_historyElementUpdate_old);
  readEigenMatrix(in, m_historyElementUpdate_new);
  readStdVector(in, m_historyElementConstPositions);
  readEigenMatrix(in, m_historyElementConst);
  readStdVector(in, m_historyMaterialUpdatePositions);
  readEigenMatrix(in, m_historyMaterialUpdate_old);
  readEigenMatrix(in, m_historyMaterialUpdate_new);
  readStdVector(in, m_historyMaterialConstPositions);
  readEigenMatrix(in, m_historyMaterialConst);
}

} // namespace HierAMuS