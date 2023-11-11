// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once


#include "shapefunctions/IntegrationsPoints/dataClasses/GaussPoints1D.h"
#include "shapefunctions/IntegrationsPoints/dataClasses/GaussPoints2D.h"
#include "shapefunctions/IntegrationsPoints/dataClasses/GaussPoints2DTriangle.h"
#include "shapefunctions/IntegrationsPoints/dataClasses/GaussPoints3D.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <vector>


namespace HierAMuS {
struct InfoData;
class IntegrationPointsManagement
{
private:    
    static GaussPoints1D GP1D;
    static GaussPoints2D GP2D;
    static GaussPoints2DTriangle GP2DTri;
    static GaussPoints3D GP3D;

    
public:
    IntegrationPointsManagement();
    ~IntegrationPointsManagement();

    static GaussPoints1D *getGaussPoints1D() { return &IntegrationPointsManagement::GP1D; };
    static GaussPoints2D *getGaussPoints2D() { return &IntegrationPointsManagement::GP2D; };
    static GaussPoints2DTriangle *getGaussPoints2DTriangle() { return &IntegrationPointsManagement::GP2DTri; };
    static GaussPoints3D *getGaussPoints3D() { return &IntegrationPointsManagement::GP3D; };
    
    static IntegrationPoints getIntegrationsPoints(indexType elementId) {
      static IntegrationPointsManagement ip;
      return IntegrationPoints(&ip, elementId);
    };
};

static IntegrationPointsManagement StaticIntegrationsPoints;

} // End Namespace