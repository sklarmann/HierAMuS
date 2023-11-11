// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <cmath>
#include "GaussPoints1D.h"
#include "GaussPoints2D.h"

namespace HierAMuS {


GaussPoints2D::GaussPoints2D(GaussPoints1D &gauss_1d) : m_gauss_1d(gauss_1d) {}


GaussPoints2D::~GaussPoints2D() {}


prec GaussPoints2D::getXi(indexType anz, indexType num) {
  indexType resNum = num % anz;
  return m_gauss_1d.getXi(anz,resNum);
}


prec GaussPoints2D::getEta(indexType anz, indexType num) {
  indexType resNum = num / anz;
  return m_gauss_1d.getXi(anz,resNum);
}


prec GaussPoints2D::getWeight(indexType anz, indexType num) {
  indexType resNum1 = num % anz;
  indexType resNum2 = num / anz;
  return m_gauss_1d.getWeight(anz,resNum1)*m_gauss_1d.getWeight(anz,resNum2);
}

} // namespace HierAMuS
