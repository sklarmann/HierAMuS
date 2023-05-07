// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

//#include "control/bindControl.hpp"
//#include "equations/bindEquation.hpp"
//#include "geometry/bindGeometry.hpp"
//#include "materials/bindMaterial.hpp"
//#include "plot/bindPlot.hpp"
//#include "solver/bindSolver.hpp"

//#include "pointercollection/PointerCollection.hpp"

#include <pybind11/pybind11.h>

#include "control/controlAddClass.h"
#include "elementFormulations/elementFormulationsAddClass.h"
#include "equations/equationsAddClass.h"
#include "finiteElements/finiteElementsAddClass.h"
#include "geometry/geometryAddClass.h"
#include "loads/loadsAddClass.h"
#include "pointercollection/PointerCollectionWrapper.h"
#include "solver/solverAddClass.h"
#include "materials/materialAddClass.h"
#include "plot/plotAddClass.h"
#include "shapefunctions/shapefunctionsAddClass.h"

namespace py = pybind11;

PYBIND11_MODULE(HierAMuSPyWrapper, m) {
  m.doc() = "pybind11 example plugin"; // optional module docstring

  auto HierAMuSPyFEMModule = m.def_submodule("HierAMuSPyFEM", "namespace");

  HierAMuS::controlAddClass control(HierAMuSPyFEMModule);
  HierAMuS::equationsAddClass equations(HierAMuSPyFEMModule);
  auto FiniteElementNameSpace =
      HierAMuSPyFEMModule.def_submodule("FiniteElement", "Finite Elements");
  HierAMuS::finiteElementsAddClass finiteElements(FiniteElementNameSpace);

  HierAMuS::PointerCollectionWrapper PointerCollection(HierAMuSPyFEMModule);
  HierAMuS::solverAdder solvers(HierAMuSPyFEMModule);
  auto GeomModule =
      HierAMuSPyFEMModule.def_submodule("Geometry", "Geometry module");
  HierAMuS::geometryAddClass geom(GeomModule);

  auto ElementFormulationNameSpace = HierAMuSPyFEMModule.def_submodule(
      "Elementformulations", "Element Formulations module");
  HierAMuS::elementFormulationsAddClass elemFormulations(ElementFormulationNameSpace);

  HierAMuS::loadsAddClass Loads(HierAMuSPyFEMModule);

  auto MaterialNameSpace =
      HierAMuSPyFEMModule.def_submodule("Materials", "Material module");
  HierAMuS::materialAddClass materials(MaterialNameSpace);

    HierAMuS::plotAddClass plot(HierAMuSPyFEMModule);

    HierAMuS::shapeFunctionsAddClass shapes(HierAMuSPyFEMModule);


  // PointerCollection.registerSelf(HierAMuSPyFEMModule);

  // PointerCollectionToPybind(HierAMuSPyFEMModule);

  // addControlFolder(HierAMuSPyFEMModule);
  // addEquationFolder(HierAMuSPyFEMModule);
  // addGeometryFolder(HierAMuSPyFEMModule);
  //addMaterialFolder(HierAMuSPyFEMModule);
  //addPlotFolder(HierAMuSPyFEMModule);
  // addSolverFolder(HierAMuSPyFEMModule);

  PointerCollection.registerFunctions();
  solvers.registerFunctions();
  geom.registerFunctions();
  equations.registerFunctions();
  control.registerFunctions();
  finiteElements.registerFunctions();
  Loads.registerFunctions();
  elemFormulations.registerFunctions();
  materials.registerFunctions();
  plot.registerFunctions();
  shapes.registerFunctions();
}