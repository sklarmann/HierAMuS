// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "forwarddeclaration.h"

#include <vector>


#include "IntegrationPointsBase.h"

namespace HierAMuS {


class GaussPoints2D : public IntegrationPointsBase
{
private:
    std::vector<std::vector<prec>> xi, eta, weight;
    indexType maxGP;
public:
    GaussPoints2D(/* args */);
    ~GaussPoints2D();

    void set1DGP(GaussPoints1D &Gauss1D);

    prec getXi(indexType anz, const indexType num);
    prec getEta(indexType anz, const indexType num);
    prec getWeight(indexType anz, indexType num);
};

}
