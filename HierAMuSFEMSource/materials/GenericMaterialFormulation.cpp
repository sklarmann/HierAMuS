// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause




#include <materials/GenericMaterialFormulation.h>
#include <control/ParameterList.h>

namespace HierAMuS::Materials {

GenericMaterialFormulation::GenericMaterialFormulation(
    PointerCollection *ptrCol) {
}

GenericMaterialFormulation::~GenericMaterialFormulation() = default;

void GenericMaterialFormulation::readData(PointerCollection& pointers, ParameterList &list) {}

auto GenericMaterialFormulation::getHistoryDataStructure(PointerCollection& pointers)
-> const HistoryDataStructure & {
  return m_historyDataStructure;
}

const HistoryDataStructure GenericMaterialFormulation::m_historyDataStructure({},{});

//HistoryDataStructure GenericMaterialFormulation::historyDataStructure =
//    HistoryDataStructure({}, {});

} // namespace HierAMuS
