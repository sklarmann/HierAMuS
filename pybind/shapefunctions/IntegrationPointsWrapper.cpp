// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "IntegrationPointsWrapper.h"
#include "pointercollection/pointercollection.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

namespace HierAMuS {

void IntegrationPointsWrapper::registerFunctions() {
  this->temp.def(py::init<IntegrationPointsManagement*,indexType>())
  .def("getTotalGP",&IntegrationPoints::getTotalGP)
  .def("getXi",&IntegrationPoints::getXi)
  .def("getEta",&IntegrationPoints::getEta)
  .def("getZeta",&IntegrationPoints::getZeta)
  .def("getWeight",&IntegrationPoints::getWeight)
  .def("setTypeOrder",&IntegrationPoints::setTypeOrder)
  .def("setCurrNumber",&IntegrationPoints::setCurrNumber)	
  .def("getIntegrationPoint",&IntegrationPoints::getIntegrationPoint,py::return_value_policy::reference)
  ;

  this->types
      .value("Gauss1D",IntegrationType::Gauss1D)
      .value("Gauss2D",IntegrationType::Gauss2D)
      .value("Gauss3D",IntegrationType::Gauss3D)
      ;
}
} // namespace HierAMuS
