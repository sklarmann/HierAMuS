// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include "IntegrationPointsManagement.h"

namespace HierAMuS {


IntegrationPointsManagement::IntegrationPointsManagement()
= default;


IntegrationPointsManagement::~IntegrationPointsManagement(/* args */)
= default;

GaussPoints1D IntegrationPointsManagement::GP1D = GaussPoints1D();
GaussPoints2D IntegrationPointsManagement::GP2D = GaussPoints2D(IntegrationPointsManagement::GP1D);
GaussPoints2DTriangle IntegrationPointsManagement::GP2DTri = GaussPoints2DTriangle();
GaussPoints3D IntegrationPointsManagement::GP3D = GaussPoints3D(IntegrationPointsManagement::GP1D);

} // End Namespace
