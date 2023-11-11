// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <elementFormulations/GenericElementFormulation.h>
#include <elementFormulations/EL101_Bernoulli2D.h>

#include <finiteElements/Edge.h>

#include <geometry/GeometryBaseData.h>
#include <geometry/VertexData.h>



#include "control/ParameterList.h"

#include <pointercollection/pointercollection.h>

#include <Eigen/Dense>

namespace HierAMuS::Elementformulations {

EL101_Bernoulli2D::EL101_Bernoulli2D(PointerCollection *ptrCol)
    : GenericElementFormulationInterface(ptrCol) {}

EL101_Bernoulli2D::~EL101_Bernoulli2D() = default;

void EL101_Bernoulli2D::readData(PointerCollection &pointers,
                                 ParameterList &list) {

  this->meshIdDisp = list.getIndexVal("meshiddisp");
  this->EA = list.getPrecVal("ea");
  this->EI = list.getPrecVal("ei");

  this->messageUnprocessed(pointers, list, "EL101_Bernoulli2D");

}

void EL101_Bernoulli2D::setDegreesOfFreedom(PointerCollection &pointers,
                                            FiniteElement::Edge &elem) {
  elem.setH1Shapes(pointers, this->meshIdDisp, 1);
}

void EL101_Bernoulli2D::AdditionalOperations(
    PointerCollection &pointers, FiniteElement::Edge &elem) {}

auto EL101_Bernoulli2D::getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
-> std::vector<DegreeOfFreedom *> {
  std::vector<DegreeOfFreedom *> Dofs;
  elem->getH1Dofs(pointers, Dofs, this->meshIdDisp, 1);
  return Dofs;
}

void EL101_Bernoulli2D::setTangentResidual(
  PointerCollection& pointers,
  FiniteElement::Edge &elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  auto &vert1 = elem.getVertex(pointers, 0);
  auto &vert2 = elem.getVertex(pointers, 1);
  Types::Vector3<prec> coor1, coor2;
  coor1 = vert1.getCoordinates();
  coor2 = vert2.getCoordinates();

  prec length;
  Dofs.clear();
  elem.getH1Dofs(pointers, Dofs, this->meshIdDisp, 1);
  IntegrationPoint ip;
  ip.xi = 0;
  length = elem.getJacobian(pointers, ip);
  length *= prec(2);

  stiffness.resize(6, 6);
  residual.resize(6);
  residual.setZero();

  stiffness.setZero();

  stiffness(0, 0) = this->EA / length;
  stiffness(0, 3) = -this->EA / length;

  stiffness(1, 1) = this->EI / length / length / length * (prec)12;
  stiffness(1, 2) = this->EI / length / length * (prec)6;
  stiffness(1, 4) = -this->EI / length / length / length * (prec)12;
  stiffness(1, 5) = this->EI / length / length * (prec)6;

  stiffness(2, 1) = this->EI / length / length * (prec)6;
  stiffness(2, 2) = this->EI / length * (prec)4;
  stiffness(2, 4) = -this->EI / length / length * (prec)6;
  stiffness(2, 5) = this->EI / length * (prec)2;

  stiffness(3, 0) = -this->EA / length;
  stiffness(3, 3) = this->EA / length;

  stiffness(4, 1) = -this->EI / length / length / length * (prec)12;
  stiffness(4, 2) = -this->EI / length / length * (prec)6;
  stiffness(4, 4) = this->EI / length / length / length * (prec)12;
  stiffness(4, 5) = -this->EI / length / length * (prec)6;

  stiffness(5, 1) = this->EI / length / length * (prec)6;
  stiffness(5, 2) = this->EI / length * (prec)2;
  stiffness(5, 4) = -this->EI / length / length * (prec)6;
  stiffness(5, 5) = this->EI / length * (prec)4;

  Eigen::Matrix<prec, 2, 1> svec;
  svec(0) = coor2[0] - coor1[0];
  svec(1) = coor2[1] - coor1[1];

  svec = svec / length;

  prec css = svec(0);
  prec sss = svec(1);

  Eigen::Matrix<prec, 6, 6> T;
  T.setZero();
  T(0, 0) = css;
  T(0, 1) = sss;
  T(1, 0) = -sss;
  T(1, 1) = css;
  T(2, 2) = 1;

  T(3, 3) = css;
  T(3, 4) = sss;
  T(4, 3) = -sss;
  T(4, 4) = css;
  T(5, 5) = 1;

  stiffness = T.transpose() * stiffness * T;
  residual = T.transpose() * residual;
}
auto EL101_Bernoulli2D::getNumberOfIntergrationPoints(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem)
    -> indexType {
  return 0;
}
void EL101_Bernoulli2D::toParaviewAdaper(
    PointerCollection &pointers, FiniteElement::Edge &elem,
    vtkPlotInterface &paraviewAdapter, ParaviewSwitch control)
{
  pointers.getSPDLogger().warn("WARNING: Plot functionality for element formulation EL101_Bernoulli2D not implemented");
}
} // End Namespace

