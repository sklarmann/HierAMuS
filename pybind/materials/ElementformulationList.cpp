// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "materials/ElementformulationList.h"
#include "pointercollection/pointercollection.h"
#include "control/ParameterList.h"

class PyElementFormulationList : public HierAMuS::Materials::ElementFormulationList {
public:
  using HierAMuS::Materials::ElementFormulationList::ElementFormulationList;

  //void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

void ElementFormulationListToPybind(py::module &m) {
  py::class_<HierAMuS::Materials::ElementFormulationList,
             //PyElementFormulationList,
             std::shared_ptr<HierAMuS::Materials::ElementFormulationList>>(
      m, "ElementFormulationList")
      .def(py::init<>())
      .def(
          "addElementFormulation",
          py::overload_cast<HierAMuS::PointerCollection&,indexType,indexType,HierAMuS::ParameterList&>(&HierAMuS::Materials::ElementFormulationList::
                                    addElementFormulation))
      
	;
}