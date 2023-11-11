// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include <elementFormulations/GenericElementFormulation.h>
#include <pointercollection/pointercollection.h>
#include <control/HandlingStructs.h>
#include "control/ParameterList.h"

#include <Eigen/Dense>

#include "finiteElements/GenericFiniteElement.h"
#include "materials/GenericMaterialFormulation.h"

#include "spdlog/fmt/ostr.h"

namespace HierAMuS::Elementformulations {

GenericElementFormulation::GenericElementFormulation(
    PointerCollection *ptrCol) {
}

GenericElementFormulation::~GenericElementFormulation() {
}





auto GenericElementFormulation::getHistoryDataStructure()
    -> const HistoryDataStructure & {
  return m_historyDataStructure;
}



void GenericElementFormulation::messageUnprocessed(PointerCollection &pointers,
                                                   ParameterList &paraMap,
                                                   std::string elementName) {

  if(!paraMap.empty()){
    auto &Logger = pointers.getSPDLogger();
    Logger.warn("nUnprocessed input parameters in elementformulation {}", elementName);

    for(auto & it : paraMap){
      Logger.warn("Parameter {} = {}", it.first, it.second);
    }
  }

}


const HistoryDataStructure GenericElementFormulation::m_historyDataStructure({},{});
}

