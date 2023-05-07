// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>


namespace py = pybind11;

#include "elementFormulations/EL290_PythonElement.h"
#include "finiteElements/GenericFiniteElement.h"
#include "equations/DegreeOfFreedom.h"
#include <vector>

namespace HierAMuS {
namespace Elementformulations {
class PyEL290_2DPythonElement : public EL290_2DPythonElement {
public:
  using EL290_2DPythonElement::EL290_2DPythonElement;

  void setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) {
    PYBIND11_OVERLOAD(void, EL290_2DPythonElement, setDegreesOfFreedom, pointers, elem);
  }

  void AdditionalOperations(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem){
    PYBIND11_OVERLOAD(void, EL290_2DPythonElement,AdditionalOperations,pointers,elem);
  }

  std::tuple<Types::MatrixXX<prec>, Types::VectorX<prec>,
             std::vector<DegreeOfFreedom *>>
  PySetTangentResidual(FiniteElement::GenericFiniteElement *elem){
    typedef std::tuple<Types::MatrixXX<prec>, Types::VectorX<prec>,
                       std::vector<DegreeOfFreedom *>> returnType;
    PYBIND11_OVERLOAD(returnType,EL290_2DPythonElement, PySetTangentResidual,elem);

  }

  std::tuple<int,int>  testfunc() {
    typedef std::tuple<int,int> returnType;
    PYBIND11_OVERLOAD(returnType ,EL290_2DPythonElement,testfunc);

  }

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class EL290_2DPythonElementWrapper {
public:
  EL290_2DPythonElementWrapper(py::module &m) : temp(m, "EL290_2DPythonElement"){};
  void registerFunctions();

private:
  typedef py::class_<EL290_2DPythonElement,PyEL290_2DPythonElement,GenericElementFormulation,
                     std::shared_ptr<EL290_2DPythonElement>>
      pw;
  pw temp;
};
} // namespace FiniteElement
} // namespace HierAMuS