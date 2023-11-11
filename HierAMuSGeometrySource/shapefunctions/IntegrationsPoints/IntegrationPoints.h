// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "datatypes.h"

#include "IntegrationPointsIterator.h"



namespace HierAMuS {

class IntegrationPointsManagement;
class IntegrationPointsBase;

enum IntegrationType { Gauss1D, Gauss2D, Gauss3D, Scaled2D, Scaled3D, Gauss2DTriangle };

struct IntegrationPoint{
  prec xi, eta, zeta, weight;
  indexType sectionNumber, gpNumber,elementId;
};

class IntegrationPoints {
private:
  IntegrationPointsManagement *IntegrationMain;
  IntegrationPointsBase *CurrIntegrationPoints{};
  indexType totGP, GpDir1, GpDir2, GpDir3;
  indexType m_currNum;

  indexType m_currSection;
  indexType m_NumberOfSections;

  IntegrationPoint point;

  indexType elementId;
  IntegrationType type;

public:
  IntegrationPoints(IntegrationPointsManagement *IntPoints, indexType elementIdIn);
  ~IntegrationPoints();

  void setType(IntegrationType type);
  void setNumGpPerDir(indexType numGP);
  void setNumGpPerDir(indexType dir1, indexType dir2);
  void setNumGpPerDir(indexType dir1, indexType dir2, indexType dir3);

  IntegrationPoint &getIntegrationPoint();

  void setTypeOrder(IntegrationType intType, indexType order);
  void setOrder(indexType order);
  void setNumberOfSections(indexType numSections);

  indexType getElementId();


  prec getXi(indexType num);
  prec getEta(indexType num);
  prec getZeta(indexType num);
  prec getWeight(indexType num);
  void setCurrNumber(indexType num);
  indexType getTotalGP();


  // Iterator
  void next();

  IntegrationPointsIterator begin();
  IntegrationPointsIterator end();


};


} // namespace HierAMuS
