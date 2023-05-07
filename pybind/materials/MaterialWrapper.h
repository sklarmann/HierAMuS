// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "materials/Material.h"

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

class MaterialWrapper {
public:
  MaterialWrapper(py::module &m)
      : temp(m, "Material"){};
  void registerFunctions();

private:
  typedef py::class_<Material, 
                     std::shared_ptr<Material>>
      pw;
  pw temp;
};
} // namespace Materials
} // namespace HierAMuS