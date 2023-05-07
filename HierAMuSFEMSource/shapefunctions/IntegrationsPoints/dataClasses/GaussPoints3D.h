// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "forwarddeclaration.h"

#include <vector>


#include "IntegrationPointsBase.h"

namespace HierAMuS {


class GaussPoints3D : public IntegrationPointsBase
{
private:
    std::vector<std::vector<prec>> xi, eta, zeta, weight;
    indexType maxGP;
public:
    GaussPoints3D(/* args */);
    ~GaussPoints3D();

    void set1DGP(GaussPoints1D &Gauss1D);

    prec getXi(indexType anz, indexType num) override;
    prec getEta(indexType anz, indexType num) override;
    prec getZeta(indexType anz, indexType num) override;
    prec getWeight(indexType anz, indexType num) override;
};

}
