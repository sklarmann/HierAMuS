// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "forwarddeclaration.h"


#include "IntegrationPointsManagement.h"

namespace HierAMuS {


IntegrationPointsManagement::IntegrationPointsManagement()
= default;


IntegrationPointsManagement::~IntegrationPointsManagement(/* args */)
= default;

void IntegrationPointsManagement::readData(InfoData &data) {
  //this->GP1D.readData(data);
  HierAMuS::IntegrationPointsManagement::GP2D.set1DGP(HierAMuS::IntegrationPointsManagement::GP1D);
  HierAMuS::IntegrationPointsManagement::GP3D.set1DGP(HierAMuS::IntegrationPointsManagement::GP1D);
}

GaussPoints1D IntegrationPointsManagement::GP1D = GaussPoints1D();
GaussPoints2D IntegrationPointsManagement::GP2D = GaussPoints2D();
GaussPoints2DTriangle IntegrationPointsManagement::GP2DTri = GaussPoints2DTriangle();
GaussPoints3D IntegrationPointsManagement::GP3D = GaussPoints3D();

} // End Namespace
