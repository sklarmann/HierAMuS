// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "solver/Constraints/BaseConstraint.h"
#include "solver/Constraints/GeneralLink.h"
#include "solver/Constraints/ConstraintHandler.h"

namespace HierAMuS {
// class PyEquationHandler : public EquationHandler {
// public:
//   using GenericFiniteElement::GenericFiniteElement;
//
//   // void renew() {
//   //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//   //}
// };

class ConstraintWrapper {
public:
  ConstraintWrapper(py::module &m)
      : ConstraintTypesModule(m, "ConstraintTypes"), BaseConstraintModule(m, "BaseConstraint"),
        GeneralLinkModule(m, "GeneralLink"), ConstraintHandlerModule(m, "ConstraintHandler"){};
  void registerFunctions();

private:
  typedef py::enum_<HierAMuS::ConstraintTypes> ConstraintWarperType;
  ConstraintWarperType ConstraintTypesModule;
  typedef py::class_<BaseConstraint, std::shared_ptr<BaseConstraint>>
      BaseConstraintType;
  BaseConstraintType BaseConstraintModule;
  typedef py::class_<GeneralLink, std::shared_ptr<GeneralLink>> GeneralLinkType;
  GeneralLinkType GeneralLinkModule;
  typedef py::class_<ConstraintHandler, std::shared_ptr<ConstraintHandler>>
      ConstraintHandlerType;
  ConstraintHandlerType ConstraintHandlerModule;
};
} // namespace HierAMuS