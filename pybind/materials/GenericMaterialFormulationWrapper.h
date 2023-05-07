// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "materials/GenericMaterialFormulation.h"

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

class GenericMaterialFormulationWrapper {
public:
  GenericMaterialFormulationWrapper(py::module &m)
      : temp(m, "GenericMaterialFormulation"), transdata(m,"MaterialTransferData"){};
  void registerFunctions();

private:
  typedef py::class_<GenericMaterialFormulation,
                     std::shared_ptr<GenericMaterialFormulation>>
      pw;
  pw temp;
  typedef py::class_<MaterialTransferData,std::shared_ptr<MaterialTransferData>> pw2;
  pw2 transdata;
};
} // namespace Materials
} // namespace HierAMuS