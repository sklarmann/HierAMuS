// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once
#include "datatypes.h"

#include "types/MatrixTypes.h"

#include "solver/HistoryDataNew/HistoryDataStructure.h"

namespace HierAMuS {
struct IntegrationPoint;
class PointerCollection;
class HistoryDataIterator;
class ParameterList;
namespace Materials {


struct MaterialTransferData {
  Types::VectorX<prec> strains;
  Types::VectorX<prec> stresses;
  Types::MatrixXX<prec> materialTangent;
  HistoryDataIterator *historyData;
};


class GenericMaterialFormulation {
public:
  explicit GenericMaterialFormulation(PointerCollection *ptrCol);
  virtual ~GenericMaterialFormulation();

  virtual void readData(PointerCollection& pointers, ParameterList &list);
  virtual void getMaterialData(PointerCollection& pointers, MaterialTransferData &inoutData, IntegrationPoint& ip) {
    HierAMuS::Materials::GenericMaterialFormulation::throwError("Error when calling Material method getMaterialData, not "
                     "implemented for given Material!");
  };

  
  virtual auto getInternalVariables(PointerCollection& pointers, MaterialTransferData &inoutData) -> std::map<std::string,Types::VectorX<prec>>{return {};};

  virtual auto getHistoryDataStructure(PointerCollection& pointers) -> const HistoryDataStructure &;

  // RVE functionality
  virtual void initRVE(PointerCollection &pointers, IntegrationPoint &ip){};
  virtual void setRVE(PointerCollection &pointers, PointerCollection &RVE) {};
  virtual void updateRVEHistory(PointerCollection &pointers, IntegrationPoint& ip) {};

protected:
  //PointerCollection *ptrCol;

private:
  static void throwError(const std::string &msg) {
    throw std::runtime_error(msg);
  }

  const static HistoryDataStructure m_historyDataStructure;
};

}
}// namespace HierAMuS
