// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <cmath>
#include "GaussPoints1D.h"
#include "GaussPoints3D.h"

namespace HierAMuS {


GaussPoints3D::GaussPoints3D(/* args */) {}


GaussPoints3D::~GaussPoints3D() {}


void GaussPoints3D::set1DGP(GaussPoints1D &Gauss1D) {
  this->maxGP = Gauss1D.getMaxGP();
  this->xi.resize(maxGP);
  this->eta.resize(maxGP);
  this->zeta.resize(maxGP);
  this->weight.resize(maxGP);

  for (auto i = 0; i < this->maxGP; ++i) {
    indexType npgs = (i + 1) * (i + 1) * (i + 1);
    this->xi[i].resize(npgs);
    this->eta[i].resize(npgs);
    this->zeta[i].resize(npgs);
    this->weight[i].resize(npgs);
    indexType pos = 0;
    for (auto n1 = 0; n1 <= i; ++n1) {
      for (auto n2 = 0; n2 <= i; ++n2) {
        for (auto n3 = 0; n3 <= i; ++n3) {
          this->xi[i][pos] = Gauss1D.getXi(i + 1, n3);
          this->eta[i][pos] = Gauss1D.getXi(i + 1, n2);
          this->zeta[i][pos] = Gauss1D.getXi(i + 1, n1);
          this->weight[i][pos] = Gauss1D.getWeight(i + 1, n1) *
                                 Gauss1D.getWeight(i + 1, n2) *
                                 Gauss1D.getWeight(i + 1, n3);
          ++pos;
        }
      }
    }
  }
}


prec GaussPoints3D::getXi(indexType anz, indexType num) {
  return this->xi[anz - 1][num];
}


prec GaussPoints3D::getEta(indexType anz, indexType num) {
  return this->eta[anz - 1][num];
}


prec GaussPoints3D::getZeta(indexType anz, indexType num) {
  return this->zeta[anz - 1][num];
}


prec GaussPoints3D::getWeight(indexType anz, indexType num) {
  return this->weight[anz - 1][num];
}

} // namespace HierAMuS
