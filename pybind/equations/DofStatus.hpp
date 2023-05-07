// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "equations/DofStatus.h"

void DofStatusToPybind(pybind11::module &m) {
  py::enum_<HierAMuS::dofStatus>(m, "dofStatus")
      .value("active", HierAMuS::dofStatus::active)
      .value("inactive", HierAMuS::dofStatus::inactive)
      .value("linked", HierAMuS::dofStatus::linked);
}
