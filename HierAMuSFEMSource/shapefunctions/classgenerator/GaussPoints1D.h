// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "forwarddeclaration.h"
#include "IntegrationPointsBase.h"

#include <vector>


namespace HierAMuS {

 
class GaussPoints1D : public IntegrationPointsBase {
private:
  static std::vector<std::vector<prec>> xi;
  static std::vector<std::vector<prec>> wi;
  indexType maxGP = 0;
  bool init;

public:
  GaussPoints1D();
  ~GaussPoints1D();

  void readData(InfoData &data);

  prec getXi(indexType anz, indexType num);
  prec getWeight(indexType anz, indexType num);

  indexType getMaxGP() { return this->xi.size(); };
};

} // namespace HierAMuS
