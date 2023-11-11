// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>
#include "IntegrationPointsBase.h"

namespace HierAMuS {

class GaussPoints2DTriangle : public IntegrationPointsBase
{
private:
    static std::vector<std::vector<prec>> xi, eta, weight;
public:
  GaussPoints2DTriangle();
  ~GaussPoints2DTriangle();

  prec getXi(indexType anz, indexType num);
  prec getEta(indexType anz, indexType num);
  prec getWeight(indexType anz, indexType num);

  indexType getMaxGP() { return static_cast<indexType>(this->xi.size()); };
};
}