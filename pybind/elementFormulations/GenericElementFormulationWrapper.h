// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "elementFormulations/GenericElementFormulation.h"

namespace HierAMuS {
namespace Elementformulations {
class PyGenericElementFormulation : public GenericElementFormulation {
public:
  using GenericElementFormulation::GenericElementFormulation;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class GenericElementFormulationWrapper {
public:
  GenericElementFormulationWrapper(py::module &m) : temp(m, "GenericElementFormulation"){};
  void registerFunctions();

private:
  typedef py::class_<GenericElementFormulation,PyGenericElementFormulation,
                     std::shared_ptr<GenericElementFormulation>>
      pw;
  pw temp;
};
} // namespace FiniteElement
} // namespace HierAMuS