// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>
#include "IntegrationPointsBase.h"

namespace HierAMuS {


class GaussPoints2D : public IntegrationPointsBase
{
private:
    GaussPoints1D &m_gauss_1d;

  public:
    GaussPoints2D(GaussPoints1D &gauss_1d);
    ~GaussPoints2D();

    prec getXi(indexType anz, const indexType num) override;
    prec getEta(indexType anz, const indexType num) override;
    prec getWeight(indexType anz, indexType num) override;
};

}
