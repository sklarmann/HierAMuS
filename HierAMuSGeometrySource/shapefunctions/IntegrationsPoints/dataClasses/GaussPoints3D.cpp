// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <cmath>
#include "GaussPoints1D.h"
#include "GaussPoints3D.h"

namespace HierAMuS {


GaussPoints3D::GaussPoints3D(GaussPoints1D &gauss_1d) : m_gauss_1d(gauss_1d) {}


GaussPoints3D::~GaussPoints3D() {}


prec GaussPoints3D::getXi(indexType anz, indexType num) {
  return m_gauss_1d.getXi(anz, num % anz);
}


prec GaussPoints3D::getEta(indexType anz, indexType num) {
  indexType resNum = (num / anz) % anz;
  return m_gauss_1d.getXi(anz, resNum);
}


prec GaussPoints3D::getZeta(indexType anz, indexType num) {
  indexType resNum = num / (anz * anz);
  return m_gauss_1d.getXi(anz, resNum);
}


prec GaussPoints3D::getWeight(indexType anz, indexType num) {
  return m_gauss_1d.getWeight(anz, num % anz) *
      m_gauss_1d.getWeight(anz, (num / anz) % anz) *
      m_gauss_1d.getWeight(anz, num / (anz * anz));
}

} // namespace HierAMuS
