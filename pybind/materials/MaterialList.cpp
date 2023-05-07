// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "materials/MaterialList.h"
#include "pointercollection/pointercollection.h"


class PyMaterialList : public HierAMuS::Materials::MaterialList {
public:
  using HierAMuS::Materials::MaterialList::MaterialList;

  //void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

void MaterialListToPybind(py::module &m) {
  py::class_<HierAMuS::Materials::MaterialList,
			 //PyMaterialList,
             std::shared_ptr<HierAMuS::Materials::MaterialList>>(
      m, "MaterialList")
      .def(py::init<>());
}