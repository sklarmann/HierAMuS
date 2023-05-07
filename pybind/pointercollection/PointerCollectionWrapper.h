// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "pointercollection/pointercollection.h"
#include "geometry/GeometryData.h"
#include "control/HandlingStructs.h"

namespace HierAMuS {
class PyPointerCollection : public HierAMuS::PointerCollection {
public:
  using HierAMuS::PointerCollection::PointerCollection;

  void renew() {
    PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  }
  std::shared_ptr<HierAMuS::Geometry::GeometryData> getGeometryData() {
    PYBIND11_OVERRIDE(std::shared_ptr<HierAMuS::Geometry::GeometryData>,
                      HierAMuS::PointerCollection, getGeometryData);
  }
};

class PointerCollectionWrapper {
public:
  PointerCollectionWrapper(py::module &m) : temp(m, "PointerCollection"){};
  void registerFunctions();

private:
  typedef py::class_<HierAMuS::PointerCollection,
                     HierAMuS::PyPointerCollection,
                     std::shared_ptr<HierAMuS::PointerCollection>>
      pw;
  pw temp;
};
} // namespace HierAMuS