// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "datatypes.h"

namespace HierAMuS {


class IntegrationPointsBase {

private:

public:
    IntegrationPointsBase();
    ~IntegrationPointsBase();
    virtual prec getXi(indexType anz, indexType num);
    virtual prec getEta(indexType anz, indexType num);
    virtual prec getZeta(indexType anz, indexType num);
    virtual prec getWeight(indexType anz, indexType num);
    virtual indexType getMaxGP();
};


} // End Namespace

