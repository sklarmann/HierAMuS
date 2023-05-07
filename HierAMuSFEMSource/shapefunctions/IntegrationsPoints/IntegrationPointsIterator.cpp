// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include "IntegrationPointsIterator.h"
#include "IntegrationPoints.h"

namespace HierAMuS{
IntegrationPointsIterator::IntegrationPointsIterator(IntegrationPoints *ptr, indexType totGP) : m_ptr(ptr), m_totGP(totGP), m_currGP(0){
  this->m_ptr->setCurrNumber(0);
}
IntegrationPointsIterator::IntegrationPointsIterator(IntegrationPoints *ptr, indexType totGP, indexType currGP) : m_ptr(ptr), m_totGP(totGP), m_currGP(currGP) {
  this->m_ptr->setCurrNumber(0);
}

indexType IntegrationPointsIterator::getElementId() { return this->m_ptr->getElementId(); }

prec IntegrationPointsIterator::Xi() {
  return this->m_ptr->getXi(this->m_currGP);
}

IntegrationPoint &IntegrationPointsIterator::operator*(){
  return this->m_ptr->getIntegrationPoint();
}

IntegrationPointsIterator &IntegrationPointsIterator::operator++(){
  ++this->m_currGP;
  this->m_ptr->setCurrNumber(this->m_currGP);
  return *this;
}

}


