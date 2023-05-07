// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <forwarddeclaration.h>

#include <map>
#include <vector>
#include <datatypes.h>
#include <Base/FEMBase.h>

namespace HierAMuS {

 class loadList : public FEMBase{
public:
  explicit loadList();
  ~loadList();
  void setLoad(indexType propNum, indexType Dof, prec loadValue, bool add);
  prec getLoad(PointerCollection& pointers, indexType Dof);
  prec getLoadIncr(PointerCollection& pointers, indexType Dof);
  void computeLoads(PointerCollection &pointers,
                    std::vector<indexType> &ids, std::vector<prec> &loads,
                    std::vector<prec> &loadincs);
  void print(PointerCollection& pointers);

private:
  //PointerCollection *pointers;
  std::map<indexType, std::map<indexType, prec>>
      loads; // contains propfunction number, dofId and value
};
} // namespace HierAMuS

