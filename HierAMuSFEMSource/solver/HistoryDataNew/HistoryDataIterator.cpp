// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "HistoryDataIterator.h"

namespace HierAMuS {

HistoryDataIterator::HistoryDataIterator(
    Eigen::Map<Types::VectorX<prec>> &elementOld,
    Eigen::Map<Types::VectorX<prec>> &elementNew,
    Eigen::Map<Types::VectorX<prec>> &elementConst,
    Eigen::Map<Types::VectorX<prec>> &materialNew,
    Eigen::Map<Types::VectorX<prec>> &materialOld,
    Eigen::Map<Types::VectorX<prec>> &materialConst)
    : m_elementStructure(nullptr), m_elementOldMap(elementOld),
      m_elementNewMap(elementNew), m_elementUpdatePosition(0),
      m_elementConstMap(elementConst), m_elementConstPosition(0),
      m_materialStructure(nullptr), m_materialOldMap(materialOld),
      m_materialNewMap(materialNew), m_materialUpdatePosition(0),
      m_materialConstMap(materialConst), m_materialConstPosition(0) {}

void HistoryDataIterator::setElementStructure(
    const HistoryDataStructure &structure) {
  m_elementStructure = &structure;
}

void HistoryDataIterator::setMaterialStructure(
    const HistoryDataStructure &structure) {
  m_materialStructure = &structure;
}

void HistoryDataIterator::next() {
  m_elementConstPosition += m_elementStructure->getNumberOfConstValues();
  m_elementUpdatePosition += m_elementStructure->getNumberOfUpdateValues();

  m_materialConstPosition += m_materialStructure->getNumberOfConstValues();
  m_materialUpdatePosition += m_materialStructure->getNumberOfUpdateValues();
}

auto HistoryDataIterator::getFieldElementUpdateNew(indexType fieldNumber)
    -> Eigen::Map<Types::MatrixXX<prec>, 0> {
  return this->getFieldMap(fieldNumber, m_elementUpdatePosition,
                           m_elementStructure->getUpdateStructure(),
                           m_elementNewMap);
}

auto HistoryDataIterator::getFieldElementUpdateOld(indexType fieldNumber)
    -> Eigen::Map<Types::MatrixXX<prec>, 0> {
  return this->getFieldMap(fieldNumber, m_elementUpdatePosition,
                           m_elementStructure->getUpdateStructure(),
                           m_elementOldMap);
}

auto HistoryDataIterator::getFieldElementConst(indexType fieldNumber)
    -> Eigen::Map<Types::MatrixXX<prec>, 0> {
  return this->getFieldMap(fieldNumber, m_elementConstPosition,
                           m_elementStructure->getConstStructure(),
                           m_elementConstMap);
}

auto HistoryDataIterator::getFieldMaterialUpdateNew(indexType fieldNumber)
    -> Eigen::Map<Types::MatrixXX<prec>, 0> {
  return this->getFieldMap(fieldNumber, m_materialUpdatePosition,
                           m_materialStructure->getUpdateStructure(),
                           m_materialNewMap);
}

auto HistoryDataIterator::getFieldMaterialUpdateOld(indexType fieldNumber)
    -> Eigen::Map<Types::MatrixXX<prec>, 0> {
  return this->getFieldMap(fieldNumber, m_materialUpdatePosition,
                           m_materialStructure->getUpdateStructure(),
                           m_materialOldMap);
}

auto HistoryDataIterator::getFieldMaterialConst(indexType fieldNumber)
    -> Eigen::Map<Types::MatrixXX<prec>, 0> {
  return this->getFieldMap(fieldNumber, m_materialConstPosition,
                           m_materialStructure->getConstStructure(),
                           m_materialConstMap);
}

auto HistoryDataIterator::getFieldMap(
    indexType fieldNumber, indexType startPosition,
    const std::vector<std::pair<indexType, indexType>> &structure,
    Eigen::Map<Types::VectorX<prec>> &data)
    -> Eigen::Map<Types::MatrixXX<prec>, 0> {

  for (auto i = 0; i < fieldNumber; ++i) {
    startPosition += structure[i].first * structure[i].second;
  }
  indexType dimx = structure[fieldNumber].first;
  indexType dimy = structure[fieldNumber].second;

  // Eigen::Map<Types::MatrixXX<prec>, 0> map(data.data() + startPosition, dimx,
  //                                          dimy);
  Eigen::Map<Types::MatrixXX<prec>, 0> map(&data(startPosition), dimx, dimy);

  return map;
}

} // namespace HierAMuS