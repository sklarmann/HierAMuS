// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "datatypes.h"
#include "types/MatrixTypes.h"
#include "HistoryDataStructure.h"
#include "Eigen/Dense"

namespace HierAMuS {

class HistoryDataIterator {

public:
  HistoryDataIterator(Eigen::Map<Types::VectorX<prec>> &elementOld,
                      Eigen::Map<Types::VectorX<prec>> &elementNew, 
                      Eigen::Map<Types::VectorX<prec>> &elementConst,
                      Eigen::Map<Types::VectorX<prec>> &materialNew,
                      Eigen::Map<Types::VectorX<prec>> &materialOld,
                      Eigen::Map<Types::VectorX<prec>> &materialConst);

  void setElementStructure(const HistoryDataStructure& structure);
  void setMaterialStructure(const HistoryDataStructure& structure);

  void next();

  auto getFieldElementUpdateNew(indexType fieldNumber)
      -> Eigen::Map<Types::MatrixXX<prec>, 0>;
  auto getFieldElementUpdateOld(indexType fieldNumber)
      -> Eigen::Map<Types::MatrixXX<prec>, 0>;
  auto getFieldElementConst(indexType fieldNumber)
      -> Eigen::Map<Types::MatrixXX<prec>, 0>;

  
  auto getFieldMaterialUpdateNew(indexType fieldNumber)
      -> Eigen::Map<Types::MatrixXX<prec>, 0>;
  auto getFieldMaterialUpdateOld(indexType fieldNumber)
      -> Eigen::Map<Types::MatrixXX<prec>, 0>;
  auto getFieldMaterialConst(indexType fieldNumber)
      -> Eigen::Map<Types::MatrixXX<prec>, 0>;

private:
  auto getFieldMap(indexType fieldNumber, indexType startPosition,
                   const std::vector<std::pair<indexType, indexType>> &structure,
                   Eigen::Map<Types::VectorX<prec>> &data)
      -> Eigen::Map<Types::MatrixXX<prec>, 0>;

  //Eigen::VectorX<prec> &m_elementOld, &m_elementNew;
  const HistoryDataStructure *m_elementStructure;
  Eigen::Map<Types::VectorX<prec>> m_elementOldMap;
  Eigen::Map<Types::VectorX<prec>> m_elementNewMap;
  indexType m_elementUpdatePosition;

  Eigen::Map<Types::VectorX<prec>> m_elementConstMap;
  indexType m_elementConstPosition;

  const HistoryDataStructure *m_materialStructure;
  Eigen::Map<Types::VectorX<prec>> m_materialOldMap;
  Eigen::Map<Types::VectorX<prec>> m_materialNewMap;
  indexType m_materialUpdatePosition;

  Eigen::Map<Types::VectorX<prec>> m_materialConstMap;
  indexType m_materialConstPosition;
};


} // namespace HierAMuS