// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "GenericFiniteElementWrapper.h"
#include "MatrixTypes.h"
#include "finiteElements/GenericFiniteElement.h"
#include "materials/GenericMaterialFormulation.h"
#include "equations/DegreeOfFreedom.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <vector>

namespace HierAMuS {
namespace FiniteElement {

void GenericFiniteElementWrapper::registerFunctions() {
  this->temp
	.def(py::init<>())
      .def("getElementType", &GenericFiniteElement::getElementType)
      .def("setAllNodeBoundaryConditionMeshId",
           &GenericFiniteElement::setAllNodeBoundaryConditionMeshId)
      .def("setId", &GenericFiniteElement::setId)
      .def("getId", &GenericFiniteElement::getId)
      .def("setMatrial", &GenericFiniteElement::setMatrial)
      .def("getMaterial", &GenericFiniteElement::getMaterial)
      //.def("getMaterialFormulation",
      //     py::overload_cast<IntegrationPoint &>(&GenericFiniteElement::getMaterialFormulation))
      .def("getMaterialId", &GenericFiniteElement::getMaterialId)
      .def("insertStiffnessResidual",
           &GenericFiniteElement::insertStiffnessResidual)
      .def("GenericSetDegreesOfFreedom",
           &GenericFiniteElement::GenericSetDegreesOfFreedom)
      .def("GenericAdditionalOperations",
           &GenericFiniteElement::GenericAdditionalOperations)
      .def("GenericSetTangentResidual",
           &GenericFiniteElement::GenericSetTangentResidual)
      .def("GenericSetMass", &GenericFiniteElement::GenericSetMass)
      .def("setVerts", &GenericFiniteElement::setVerts)
      .def("setEdges", &GenericFiniteElement::setEdges)
      .def("setFace", &GenericFiniteElement::setFace)
      .def("setVolume", &GenericFiniteElement::setVolume)
      .def("setSpecial", py::overload_cast<std::vector<indexType>&>(&GenericFiniteElement::setSpecial))
      .def("setSpecial", py::overload_cast<indexType>(&GenericFiniteElement::setSpecial))
      .def("setH1Shapes",&GenericFiniteElement::setH1Shapes)

      .def("getH1Dofs",[](GenericFiniteElement *self,PointerCollection &pointers,indexType meshId, indexType order){
        std::vector<DegreeOfFreedom*> Dofs;
        self->getH1Dofs(pointers,Dofs,meshId,order);
        return Dofs;
      },py::return_value_policy::reference)
      .def("setAllNodeBoundaryConditionMeshId",&GenericFiniteElement::setAllNodeBoundaryConditionMeshId)

      .def("getJacobian",([](GenericFiniteElement &self,PointerCollection &pointers,prec xsi){
        prec jaco;
        self.getJacobian(pointers,jaco,xsi);
        return jaco;
      }))

      .def("getJacobian",([](GenericFiniteElement &self,PointerCollection &pointers,prec xsi,prec eta){

        IntegrationPoint ip;
        ip.xi = xsi;
        ip.eta = eta;
        Types::MatrixXX<prec> jacobi = self.getJacobian(pointers,ip);
        return jacobi;
      }))

      .def("getH1Shapes",[](GenericFiniteElement &self, PointerCollection &pointers, indexType order, Types::MatrixXX<prec> &jacobi, prec xsi, prec eta){
        Types::VectorX<prec> shape;
        Types::Matrix2X<prec> shapederiv;
        IntegrationPoint ip;
        ip.xi = xsi;
        ip.eta = eta;
        auto H1Shapes = self.getH1Shapes(pointers,order,jacobi,ip);
        return std::make_tuple(H1Shapes.shapes, H1Shapes.shapeDeriv);
      })
      .def("getMaterialFormulation",py::overload_cast<PointerCollection&>(&GenericFiniteElement::getMaterialFormulation))
      .def("getMaterialFormulation",
           py::overload_cast<PointerCollection & , IntegrationPoint &>(
               &GenericFiniteElement::getMaterialFormulation))

      .def("getSolution",[](GenericFiniteElement &self,PointerCollection &pointers,std::vector<DegreeOfFreedom*> &Dofs){
        Types::VectorX<prec> sol;
        self.getSolution(pointers,Dofs,sol);
        return sol;
      })
      .def("computeNorm",&GenericFiniteElement::computeNorm)
      .def("setFaces",&GenericFiniteElement::setFaces)
      .def("setMaterialPerSubElement",&GenericFiniteElement::setMaterialPerSubElement);



	;
}
} // namespace FiniteElement
} // namespace HierAMuS
