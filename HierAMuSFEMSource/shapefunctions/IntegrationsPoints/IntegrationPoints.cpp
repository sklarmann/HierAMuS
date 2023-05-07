// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "IntegrationPoints.h"
#include "dataClasses/GaussPoints1D.h"
#include "helperClasses/IntegrationPointsManagement.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPointsIterator.h"
#include <stdexcept>

namespace HierAMuS {

IntegrationPoints::IntegrationPoints(IntegrationPointsManagement *IntPoints, indexType elementIdIn)
    : IntegrationMain(IntPoints), m_currNum(0), totGP(0), elementId(elementIdIn){};

IntegrationPoints::~IntegrationPoints() = default;

void IntegrationPoints::setTypeOrder(IntegrationType intType, indexType order) {
  throw std::runtime_error("Deprecated Gauss Point routine called");
  // this->type = intType;

  // switch (this->type) {
  // case Gauss1D:
  //   this->CurrIntegrationPoints =
  //   HierAMuS::IntegrationPointsManagement::getGaussPoints1D(); this->GpDir1
  //   = order + 1; if (this->GpDir1 % 2 == 0) {
  //     this->GpDir1 /= 2;
  //   } else {
  //     this->GpDir1 += 1;
  //     this->GpDir1 /= 2;
  //   }
  //   this->totGP = this->GpDir1;
  //   break;
  // case Gauss2D:
  //   this->CurrIntegrationPoints =
  //   HierAMuS::IntegrationPointsManagement::getGaussPoints2D(); this->GpDir1
  //   = order + 1; if (this->GpDir1 % 2 == 0) {
  //     this->GpDir1 /= 2;
  //   } else {
  //     this->GpDir1 += 1;
  //     this->GpDir1 /= 2;
  //   }
  //   this->GpDir2 = this->GpDir1;
  //   this->totGP = this->GpDir1 * this->GpDir2;
  //   break;
  // case Gauss3D:
  //   this->CurrIntegrationPoints =
  //   HierAMuS::IntegrationPointsManagement::getGaussPoints3D(); this->GpDir1
  //   = order + 1; if (this->GpDir1 % 2 == 0) {
  //     this->GpDir1 /= 2;
  //   } else {
  //     this->GpDir1 += 1;
  //     this->GpDir1 /= 2;
  //   }
  //   this->GpDir2 = this->GpDir1;
  //   this->GpDir3 = this->GpDir1;
  //   this->totGP = this->GpDir1 * this->GpDir1 * this->GpDir1;
  //   break;
  // }
}

void IntegrationPoints::setOrder(indexType order) {
  switch (this->type) {
  case Gauss1D:
    this->GpDir1 = order + 1;
    if (this->GpDir1 % 2 == 0) {
      this->GpDir1 /= 2;
    } else {
      this->GpDir1 += 1;
      this->GpDir1 /= 2;
    }
    this->totGP = this->GpDir1;
    break;
  case Gauss2D:
    this->GpDir1 = order + 1;
    if (this->GpDir1 % 2 == 0) {
      this->GpDir1 /= 2;
    } else {
      this->GpDir1 += 1;
      this->GpDir1 /= 2;
    }
    this->GpDir2 = this->GpDir1;
    this->totGP = this->GpDir1 * this->GpDir2;
    break;
  case Gauss3D:
    this->GpDir1 = order + 1;
    if (this->GpDir1 % 2 == 0) {
      this->GpDir1 /= 2;
    } else {
      this->GpDir1 += 1;
      this->GpDir1 /= 2;
    }
    this->GpDir2 = this->GpDir1;
    this->GpDir3 = this->GpDir1;
    this->totGP = this->GpDir1 * this->GpDir1 * this->GpDir1;
    break;
  case Scaled2D:
    this->GpDir1 = order + 1;
    if (this->GpDir1 % 2 == 0) {
      this->GpDir1 /= 2;
    } else {
      this->GpDir1 += 1;
      this->GpDir1 /= 2;
    }
    this->GpDir2 = this->GpDir1;
    this->totGP = this->GpDir1 * this->GpDir2 * this->m_NumberOfSections;
    break;
  case Scaled3D:
    this->GpDir1 = order + 1;
    if (this->GpDir1 % 2 == 0) {
      this->GpDir1 /= 2;
    } else {
      this->GpDir1 += 1;
      this->GpDir1 /= 2;
    }
    this->GpDir2 = this->GpDir1;
    this->GpDir3 = this->GpDir1;
    this->totGP =
        this->GpDir1 * this->GpDir2 * this->GpDir3 * this->m_NumberOfSections;
    break;
  case Gauss2DTriangle:
    switch (order) {
      case 1:
        this->totGP = 1;
      break;
      case 2:
        this->totGP = 3;
      break;
      case 3:
        this->totGP = 4;
      break;
      case 4:
        this->totGP = 6;
      break;
      case 5:
        this->totGP = 7;
      break;
      case 6:
        this->totGP = 12;
      break;
      case 7:
        this->totGP = 13;
      break;
      case 8:
        this->totGP = 16;
      break;
      case 9:
        this->totGP = 19;
      break;
      case 10:
        this->totGP = 25;
      break;
      case 11:
        this->totGP = 27;
      break;
      case 12:
        this->totGP = 33;
      break;
      case 13:
        this->totGP = 37;
      break;
      case 14:
        this->totGP = 42;
      break;
      case 15:
        this->totGP = 48;
      break;
      case 16:
        this->totGP = 52;
      break;
      case 17:
        this->totGP = 61;
      break;
    }
    this->GpDir1 = order;
  }
}

void IntegrationPoints::setNumberOfSections(indexType numSections) {
  this->m_NumberOfSections = numSections;
}

indexType IntegrationPoints::getElementId() { return elementId; }

prec IntegrationPoints::getXi(indexType num) {

  return this->CurrIntegrationPoints->getXi(this->GpDir1, num);
}

prec IntegrationPoints::getEta(indexType num) {

  return this->CurrIntegrationPoints->getEta(this->GpDir1, num);
}

prec IntegrationPoints::getZeta(indexType num) {

  return this->CurrIntegrationPoints->getZeta(this->GpDir1, num);
}

prec IntegrationPoints::getWeight(indexType num) {

  return this->CurrIntegrationPoints->getWeight(this->GpDir1, num);
}

indexType IntegrationPoints::getTotalGP() { return this->totGP; }

void IntegrationPoints::next() { this->setCurrNumber(this->m_currNum + 1); }

IntegrationPointsIterator IntegrationPoints::begin() {
  return {this, this->totGP};
}

IntegrationPointsIterator IntegrationPoints::end() {
  return {this, this->totGP, this->totGP};
}

void IntegrationPoints::setType(IntegrationType type) {

  this->type = type;
  switch (type) {
  case IntegrationType::Gauss1D:
    this->CurrIntegrationPoints =
        HierAMuS::IntegrationPointsManagement::getGaussPoints1D();
    break;
  case IntegrationType::Gauss2D:
    this->CurrIntegrationPoints =
        HierAMuS::IntegrationPointsManagement::getGaussPoints2D();
    break;
  case IntegrationType::Gauss3D:
    this->CurrIntegrationPoints =
        HierAMuS::IntegrationPointsManagement::getGaussPoints3D();
    break;
  case IntegrationType::Scaled2D:
    this->CurrIntegrationPoints =
        HierAMuS::IntegrationPointsManagement::getGaussPoints2D();
    break;
  case IntegrationType::Scaled3D:
    this->CurrIntegrationPoints =
        HierAMuS::IntegrationPointsManagement::getGaussPoints3D();
    break;
  case IntegrationType::Gauss2DTriangle:
    this->CurrIntegrationPoints =
        HierAMuS::IntegrationPointsManagement::getGaussPoints2DTriangle();
  }
}

void IntegrationPoints::setNumGpPerDir(indexType numGP) {
  this->GpDir1 = numGP;
  switch (this->type) {
  case IntegrationType::Gauss1D:
    this->totGP = numGP;
    break;
  case IntegrationType::Gauss2D:
    this->totGP = numGP * numGP;
    break;
  case IntegrationType::Gauss3D:
    this->totGP = numGP * numGP * numGP;
    break;
  case IntegrationType::Scaled2D:
    this->totGP = numGP * numGP * this->m_NumberOfSections;
    break;
    case IntegrationType::Scaled3D:
    
    break;
  case IntegrationType::Gauss2DTriangle:
    this->totGP = numGP;
    break;
  }
}
void IntegrationPoints::setNumGpPerDir(indexType dir1, indexType dir2) {
  if (this->type == IntegrationType::Gauss2D) {
    this->GpDir1 = dir1;
    this->GpDir2 = dir2;
    this->totGP = dir1 * dir2;
  }
}

void IntegrationPoints::setCurrNumber(indexType num) {
  this->m_currNum = num;
  switch (this->type) {
  case IntegrationType::Scaled2D:
  case IntegrationType::Scaled3D: {
    indexType GPperSection = this->totGP / this->m_NumberOfSections;
    indexType section = num / GPperSection;
    indexType gp = num % GPperSection;
    this->point = {this->CurrIntegrationPoints->getXi(this->GpDir1, gp),
                   this->CurrIntegrationPoints->getEta(this->GpDir1, gp),
                   this->CurrIntegrationPoints->getZeta(this->GpDir1, gp),
                   this->CurrIntegrationPoints->getWeight(this->GpDir1, gp),
                   section,
                   m_currNum,elementId};
    break;
  }
  default: {
    if (this->m_currNum < this->totGP) {
      this->point = {
          this->CurrIntegrationPoints->getXi(this->GpDir1, m_currNum),
          this->CurrIntegrationPoints->getEta(this->GpDir1, m_currNum),
          this->CurrIntegrationPoints->getZeta(this->GpDir1, m_currNum),
          this->CurrIntegrationPoints->getWeight(this->GpDir1, m_currNum),
          0,
          m_currNum,elementId};
    }
  } break;
  }
}

} // namespace HierAMuS
