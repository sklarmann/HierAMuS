// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once
#include "datatypes.h"

namespace HierAMuS {
class IntegrationPoints;
class IntegrationPoint;
class IntegrationPointsIterator {
public:
  IntegrationPointsIterator(IntegrationPoints *ptr, indexType totGP);
  IntegrationPointsIterator(IntegrationPoints *ptr, indexType totGP, indexType currGP);

  IntegrationPointsIterator &operator++();

  IntegrationPoint & operator*();

  friend bool operator==(const IntegrationPointsIterator& a, const IntegrationPointsIterator& b){return a.m_currGP == b.m_totGP;};
  friend bool operator!=(const IntegrationPointsIterator& a, const IntegrationPointsIterator& b){return a.m_currGP != b.m_totGP;};

  indexType getElementId();

  prec Xi();

private:
  IntegrationPoints *m_ptr;
  indexType m_totGP;
  indexType m_currGP;
};
}