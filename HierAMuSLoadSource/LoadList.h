// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "datatypes.h"

#include <map>
#include <vector>


#include "spdlog/spdlog.h"

namespace HierAMuS {
class PropfunctionHandler;
class LoadList {
public:
  explicit LoadList();
  ~LoadList();
  void setLoad(indexType propNum, indexType Dof, prec loadValue, bool add);
  prec getLoad(PropfunctionHandler &propHandler, indexType Dof);
  prec getLoadIncr(PropfunctionHandler &propHandler, indexType Dof);
  void computeLoads(PropfunctionHandler &propHandler,
                    std::vector<indexType> &ids, std::vector<prec> &loads,
                    std::vector<prec> &loadincs);
  void print(spdlog::logger &pointers);

private:
  //PointerCollection *pointers;
  std::map<indexType, std::map<indexType, prec>>
      loads; // contains propfunction number, dofId and value
};
} // namespace HierAMuS

