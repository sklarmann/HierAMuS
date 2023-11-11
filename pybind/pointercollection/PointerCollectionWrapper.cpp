// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "PointerCollectionWrapper.h"

#include "EquationHandler.h"
#include "finiteElements/ElementList.h"
#include "materials/MaterialList.h"
#include "materials/ElementformulationList.h"
#include "materials/MaterialformulationList.h"
#include "LoadList.h"
#include "PropfunctionHandler.h"
#include "solver/GenericSolutionState.h"
#include "plot/plotControl.h"
#include "plot/vtkplotClass.h"
#include "control/ParameterList.h"

void HierAMuS::PointerCollectionWrapper::registerFunctions() {
  this->temp.def(py::init<>())
      .def("renew", &HierAMuS::PointerCollection::renew,
           "Resets the data in class PointerCollection.")
      // Geometry
      .def("setGeometryData", &HierAMuS::PointerCollection::setGeometryData)
      .def("getGeometryData", &HierAMuS::PointerCollection::getGeometryData,
           "Returns the class GeometryData.")
      .def("newGeometry", &HierAMuS::PointerCollection::newGeometry,
           "Deletes old GeometryData and creates a new one.")
      // Solution State
      .def("setSolutionState",
           py::overload_cast<std::shared_ptr<HierAMuS::GenericSolutionState>>(
               &HierAMuS::PointerCollection::setSolutionState))
      .def("setSolutionState",
           py::overload_cast<HierAMuS::SolutionTypes,
                             HierAMuS::ParameterList &>(
               &HierAMuS::PointerCollection::setSolutionState))
      .def("getSolutionState", &HierAMuS::PointerCollection::getSolutionState)
      // EquationHandler
      .def("setEquationHandler",
           &HierAMuS::PointerCollection::setEquationHandler)
      .def("getEquationHandler",
           &HierAMuS::PointerCollection::getEquationHandler)
      // InfoData
      .def("setInfoData", &HierAMuS::PointerCollection::setInfoData)
      // ElementsList
      .def("getElementList", &HierAMuS::PointerCollection::getElementList)
      // MaterialsList
      .def("getMaterialList", &HierAMuS::PointerCollection::getMaterialList)
      // ElementFormulationList
      .def("getElementFormulationList", &HierAMuS::PointerCollection::getElementFormulationList)
      // MaterialFormulationList
      .def("getMaterialFormulationList",&PointerCollection::getMaterialFormulationList)
    // Proploads
      .def("getPropLoads",&PointerCollection::getPropLoads)
    // vtkPlot
      .def("getVtkPlotInterface",&PointerCollection::getVtkPlotInterface)
      .def("getLoadList",&HierAMuS::PointerCollection::getLoadList)
      .def("getPrescribedDisplacements",&HierAMuS::PointerCollection::getPrescribedDisplacements)
      .def("setMaxThreads", &PointerCollection::setMaxThreads)
      .def("getPlotControlInterface",
           &PointerCollection::getPlotControlInterface)
      .def("solutionStateToFile", &PointerCollection::solutionStateToFile)
      .def("solutionStateFromFile", &PointerCollection::solutionStateFromFile)
	;
}
