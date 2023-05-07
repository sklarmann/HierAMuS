// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "forwarddeclaration.h"

#include "../dataClasses/GaussPoints1D.h"
#include "../dataClasses/GaussPoints2D.h"
#include "../dataClasses/GaussPoints2DTriangle.h"
#include "../dataClasses/GaussPoints3D.h"
#include <vector>


namespace HierAMuS {

class IntegrationPointsManagement
{
private:
    void computeGP1D();
    
    static GaussPoints1D GP1D;
    static GaussPoints2D GP2D;
    static GaussPoints2DTriangle GP2DTri;
    static GaussPoints3D GP3D;

    
public:
    IntegrationPointsManagement();
    ~IntegrationPointsManagement();
    void readData(InfoData &data);

    static GaussPoints1D *getGaussPoints1D() { return &IntegrationPointsManagement::GP1D; };
    static GaussPoints2D *getGaussPoints2D() { return &IntegrationPointsManagement::GP2D; };
    static GaussPoints2DTriangle *getGaussPoints2DTriangle() { return &IntegrationPointsManagement::GP2DTri; };
    static GaussPoints3D *getGaussPoints3D() { return &IntegrationPointsManagement::GP3D; };
    
};

} // End Namespace