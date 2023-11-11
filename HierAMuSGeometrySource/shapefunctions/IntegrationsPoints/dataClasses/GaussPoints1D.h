// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "IntegrationPointsBase.h"
#include <vector>


namespace HierAMuS {

 
class GaussPoints1D : public IntegrationPointsBase {
private:
  static std::vector<std::vector<prec>> xi;
  static std::vector<std::vector<prec>> wi;

public:
  GaussPoints1D();
  ~GaussPoints1D();


  prec getXi(indexType anz, indexType num);
  prec getWeight(indexType anz, indexType num);

  indexType getMaxGP() { return static_cast<indexType>(this->xi.size()); };
};

} // namespace HierAMuS
