// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <cmath>
#include "GaussPoints1D.h"
#include "GaussPoints2D.h"

namespace HierAMuS {


GaussPoints2D::GaussPoints2D(/* args */) {}


GaussPoints2D::~GaussPoints2D() {}


void GaussPoints2D::set1DGP(GaussPoints1D &Gauss1D) {
  this->maxGP = Gauss1D.getMaxGP();
  this->xi.resize(maxGP);
  this->eta.resize(maxGP);
  this->weight.resize(maxGP);

  for (auto i = 0; i < this->maxGP; ++i) {
    indexType npgs = (i + 1) * (i + 1);
    this->xi[i].resize(npgs);
    this->eta[i].resize(npgs);
    this->weight[i].resize(npgs);
    indexType pos = 0;
    for (auto n1 = 0; n1 <= i; ++n1) {
      for (auto n2 = 0; n2 <= i; ++n2) {
        this->xi[i][pos] = Gauss1D.getXi(i + 1, n2);
        this->eta[i][pos] = Gauss1D.getXi(i + 1, n1);
        this->weight[i][pos] =
            Gauss1D.getWeight(i + 1, n1) * Gauss1D.getWeight(i + 1, n2);
        ++pos;
      }
    }
  }
}


prec GaussPoints2D::getXi(indexType anz, indexType num) {
  return this->xi[anz-1][num];
}


prec GaussPoints2D::getEta(indexType anz, indexType num) {
  return this->eta[anz-1][num];
}


prec GaussPoints2D::getWeight(indexType anz, indexType num) {
  return this->weight[anz-1][num];
}

} // namespace HierAMuS
