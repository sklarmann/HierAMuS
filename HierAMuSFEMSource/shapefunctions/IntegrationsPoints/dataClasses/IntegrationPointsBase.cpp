// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "IntegrationPointsBase.h"

namespace HierAMuS {


IntegrationPointsBase::IntegrationPointsBase() = default;;


IntegrationPointsBase::~IntegrationPointsBase() = default;;


prec IntegrationPointsBase::getXi(indexType anz, indexType num) {return 0;};



prec IntegrationPointsBase::getEta(indexType anz, indexType num) {return 0;};



prec IntegrationPointsBase::getZeta(indexType anz, indexType num) {return 0;};



prec IntegrationPointsBase::getWeight(indexType anz, indexType num) {return 0;};



indexType IntegrationPointsBase::getMaxGP() {return 0;};



}  // End Namespace

