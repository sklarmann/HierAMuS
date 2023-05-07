// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "materials/MaterialList.h"

namespace HierAMuS {
namespace Materials {

//class PyTransientSolutionNewmark : public MaterialList {
//public:
//  using MaterialList::MaterialList;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

class MaterialListWrapper {
public:
  MaterialListWrapper(py::module &m)
      : temp(m, "MaterialList"){};
  void registerFunctions();

private:
  typedef py::class_<MaterialList, 
                     std::shared_ptr<MaterialList>>
      pw;
  pw temp;
};
} // namespace Materials
} // namespace HierAMuS