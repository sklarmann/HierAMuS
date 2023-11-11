// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

namespace HierAMuS {
//class PyEquationHandler : public EquationHandler {
//public:
//  using GenericFiniteElement::GenericFiniteElement;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};



class IntegrationPointsWrapper {
public:
  IntegrationPointsWrapper(py::module &m)
      : tempip(m, "IntegrationPoint"), temp(m, "IntegrationPoints"),
        types(m, "IntegrationType"){};
  void registerFunctions();

private:
  typedef py::class_<IntegrationPoint, std::shared_ptr<IntegrationPoint>>
      pw3;
  pw3 tempip;
  typedef py::class_<IntegrationPoints, std::shared_ptr<IntegrationPoints>>
      pw;
  pw temp;

  typedef py::enum_<IntegrationType>
    pw2;
    pw2 types;
};
} // namespace HierAMuS