// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "GenericMaterialFormulationWrapper.h"
#include "materials/GenericMaterialFormulation.h"
#include "pointercollection/pointercollection.h"
#include <pybind11/eigen.h>


namespace HierAMuS {
namespace Materials {

void GenericMaterialFormulationWrapper::registerFunctions() {

  this->temp.def("getMaterialData",&GenericMaterialFormulation::getMaterialData)
            .def("setRVE",&GenericMaterialFormulation::setRVE);




  this->transdata
    .def(py::init<>())
      .def_readwrite("strains",&MaterialTransferData::strains)
      .def_readwrite("stresses",&MaterialTransferData::stresses)
      .def_readwrite("materialTangent",&MaterialTransferData::materialTangent)
  ;

}
} // namespace Materials
} // namespace HierAMuS
