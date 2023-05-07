// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "pointercollection/pointercollection.h"
#include "solver/TransientSolutionNewmark.h"

namespace HierAMuS {
class PyTransientSolutionNewmark : public HierAMuS::TransientSolutionNewmark {
public:
  using HierAMuS::TransientSolutionNewmark::TransientSolutionNewmark;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class TransientSolutionNewmarkWrapper {
public:
  TransientSolutionNewmarkWrapper(py::module &m)
      : temp(m, "TransientSolutionNewmark"){};
  void registerFunctions();

private:
  typedef py::class_<TransientSolutionNewmark, PyTransientSolutionNewmark,
                     std::shared_ptr<TransientSolutionNewmark>>
      pw;
  pw temp;
};
} // namespace HierAMuS