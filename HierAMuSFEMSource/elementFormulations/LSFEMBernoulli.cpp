// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <elementFormulations/GenericElementFormulation.h>
#include <elementFormulations/LSFEMBernoulli.h>

#include <finiteElements/GenericFiniteElement.h>

#include <geometry/Base.h>
#include <geometry/Vertex.h>


#include <equations/NodeSet.h>
#include <equations/Nodetypes.h>
#include <equations/GenericNodes.h>
#include <equations/EquationHandler.h>

#include <pointercollection/pointercollection.h>
#include <solver/GenericSolutionState.h>

#include <Eigen/Dense>

#include <loads/PropfunctionHandler.h>

namespace HierAMuS::Elementformulations {

LSFEMBernoulli::LSFEMBernoulli(PointerCollection *ptrCol)
    : GenericElementFormulation(ptrCol) {}

LSFEMBernoulli::~LSFEMBernoulli() = default;

void LSFEMBernoulli::readData(PointerCollection &pointers,
                              ParameterList &list) {

  this->meshIdDisp = list.getIndexVal("meshiddisp");
  this->EA = list.getPrecVal("ea");
  this->EI = list.getPrecVal("ei");
  this->GA = list.getPrecVal("ga");
  this->rhoA = list.getPrecVal("rhoa");
  this->propnum = list.getIndexVal("propnum");
  this->localLoad = list.getIndexVal("local");
  this->qx = list.getPrecVal("qx");
  this->qy = list.getPrecVal("qy");
  this->mz = list.getPrecVal("mz");

  this->messageUnprocessed(pointers, list, "Element LSFEMBernoulli");
}

void LSFEMBernoulli::setDegreesOfFreedom(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) {

  for (auto i = 0; i < 2; ++i) {
    auto &geoTemp = elem->getVertex(pointers, i);
    geoTemp.setNodeSet(pointers, this->meshIdDisp, 2,
                        NodeTypes::displacement);
  }
}

void LSFEMBernoulli::setTangentResidual(
  PointerCollection& pointers,
  FiniteElement::GenericFiniteElement *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  auto eqHandler = pointers.getEquationHandler();
  NodeSet *set1, *set2;
  auto &vert1 = elem->getVertex(pointers, 0);
  auto &vert2 = elem->getVertex(pointers, 1);
  set1 = vert1.getSetMeshId(pointers, this->meshIdDisp);
  set2 = vert2.getSetMeshId(pointers, this->meshIdDisp);
  std::vector<GenericNodes *> temp, Nodes;
  eqHandler->getNodes(temp, *set1);
  Nodes.push_back(temp[0]);
  Nodes.push_back(temp[1]);
  temp.clear();
  eqHandler->getNodes(temp, *set2);
  // set2->getNodes(this->ptrCol, temp);
  Nodes.push_back(temp[0]);
  Nodes.push_back(temp[1]);

  std::vector<DegreeOfFreedom *> tempDofs;
  Nodes[0]->getDegreesOfFreedom(pointers, tempDofs);
  for (auto i = 0; i < 3; ++i) {
    Dofs.push_back(tempDofs[i]);
  }
  Nodes[1]->getDegreesOfFreedom(pointers, tempDofs);
  for (auto i = 0; i < 3; ++i) {
    Dofs.push_back(tempDofs[i]);
  }
  Nodes[2]->getDegreesOfFreedom(pointers, tempDofs);
  for (auto i = 0; i < 3; ++i) {
    Dofs.push_back(tempDofs[i]);
  }
  Nodes[3]->getDegreesOfFreedom(pointers, tempDofs);
  for (auto i = 0; i < 3; ++i) {
    Dofs.push_back(tempDofs[i]);
  }

  Types::Vector3<prec> coor1, coor2;
  coor1 = vert1.getCoordinates();
  coor2 = vert2.getCoordinates();
  prec length = (coor1 - coor2).norm();

  auto sol = pointers.getSolutionState();
  Eigen::Matrix<prec, Eigen::Dynamic, 1> disp;
  disp = sol->getSolution(Dofs);

  std::vector<prec> gp, weight;

  gp.push_back((prec)-0.5);
  gp.push_back((prec)0.5);
  weight.push_back((prec)1 / sqrt((prec)3));
  weight.push_back((prec)1 / sqrt((prec)3));

  Eigen::Matrix<prec, 2, 1> svec;
  svec(0) = coor2[0] - coor1[0];
  svec(1) = coor2[1] - coor1[1];

  svec = svec / length;

  Eigen::Matrix<prec, 4, Eigen::Dynamic> Bmat;
  Eigen::Matrix<prec, 1, Eigen::Dynamic> NN;
  Bmat.resize(4, 12);
  Bmat.setZero();
  NN.resize(1, 12);
  NN.setZero();

  stiffness.resize(12, 12);
  stiffness.setZero();
  residual.resize(12);
  residual.setZero();
  Eigen::Matrix<prec, Eigen::Dynamic, 1> shapeDeriv, shape;
  shapeDeriv.resize(2);
  shape.resize(2);
  for (auto i = 0; i < gp.size(); ++i) {
    shape(0) = ((prec)1.0 - gp[i]) / (prec)2.0;
    shape(1) = ((prec)1.0 + gp[i]) / (prec)2.0;
    shapeDeriv(0) = -(prec)1.0 / length;
    shapeDeriv(1) = (prec)1.0 / length;

    for (auto j = 0; j < 2; ++j) {
      Bmat(0, j * 6) = shapeDeriv(j);
      Bmat(0, j * 6 + 1) = -shape(j);

      Bmat(1, j * 6 + 1) = this->EI * shapeDeriv(j);
      Bmat(1, j * 6 + 2) = -shape(j);

      Bmat(2, j * 6 + 2) = shapeDeriv(j);
      Bmat(2, j * 6 + 3) = -shape(j);

      Bmat(3, j * 6 + 3) = shapeDeriv(j);

      NN(0, j * 6 + 3) = (prec)1.0 / (prec)2.0;
      ;
    }

    residual += NN.transpose() * this->qy; // *length / (prec)2.0 * weight[i];
    stiffness += (Bmat.transpose() * Bmat) * length / (prec)2.0 * weight[i];
  }
}

void LSFEMBernoulli::setMass(
  PointerCollection& pointers,
  FiniteElement::GenericFiniteElement *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {
  auto eqHandler = pointers.getEquationHandler();
  NodeSet *set1, *set2;
  auto &vert1 = elem->getVertex(pointers, 0);
  auto &vert2 = elem->getVertex(pointers, 1);
  set1 = vert1.getSetMeshId(pointers, this->meshIdDisp);
  set2 = vert2.getSetMeshId(pointers, this->meshIdDisp);
  std::vector<GenericNodes *> temp, Nodes;
  eqHandler->getNodes(temp, *set1);
  Nodes.push_back(temp[0]);
  temp.clear();
  eqHandler->getNodes(temp, *set2);
  // set2->getNodes(*this->ptrCol, temp);
  Nodes.push_back(temp[0]);

  std::vector<DegreeOfFreedom *> tempDofs;
  Nodes[0]->getDegreesOfFreedom(pointers, tempDofs);
  for (auto i = 0; i < 3; ++i) {
    Dofs.push_back(tempDofs[i]);
  }
  Nodes[1]->getDegreesOfFreedom(pointers, tempDofs);
  for (auto i = 0; i < 3; ++i) {
    Dofs.push_back(tempDofs[i]);
  }

  Types::Vector3<prec> coor1, coor2;
  coor1 = vert1.getCoordinates();
  coor2 = vert2.getCoordinates();
  prec length = (coor1 - coor2).norm();

  stiffness.resize(6, 6);
  for (auto i = 0; i < 2; ++i) {
    stiffness(i, i) = this->rhoA * length / (prec)2;
    stiffness(i + 3, i + 3) = this->rhoA * length / (prec)2;
  }
  // stiffness(2, 2) = this->rhoA*length / (prec)2*(prec)1e-12;
  // stiffness(5, 5) = this->rhoA*length / (prec)2*(prec)1e-12;
}

void LSFEMBernoulli::getElementsLocalNodalReactions(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem,
    std::map<indexType, std::vector<prec>> &vReacs) {
  indexType vert1 = elem->getVertexId(pointers, 0);
  indexType vert2 = elem->getVertexId(pointers, 1);

  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> stiffness;
  Eigen::Matrix<prec, Eigen::Dynamic, 1> residual;
  std::vector<DegreeOfFreedom *> Dofs;

  this->setTangentResidual(pointers, elem, stiffness, residual, Dofs);

  auto &gvert1 = elem->getVertex(pointers, 0);
  auto &gvert2 = elem->getVertex(pointers, 1);

  Types::Vector3<prec> coor1, coor2;
  coor1 = gvert1.getCoordinates();
  coor2 = gvert2.getCoordinates();
  prec length = (coor1 - coor2).norm();

  Eigen::Matrix<prec, 2, 1> svec;
  svec(0) = coor2[0] - coor1[0];
  svec(1) = coor2[1] - coor1[1];

  svec = svec / length;

  prec css, sss;

  css = svec(0);
  sss = svec(1);

  Eigen::Matrix<prec, 6, 6> T;
  T.setZero();
  T(0, 0) = css;
  T(0, 1) = sss;
  T(1, 0) = -sss;
  T(1, 1) = css;
  T(2, 2) = (prec)1.0;

  T(3, 3) = css;
  T(3, 4) = sss;
  T(4, 3) = -sss;
  T(4, 4) = css;
  T(5, 5) = (prec)1.0;

  residual = T * residual;

  std::vector<prec> temp1(6), temp2(6);
  for (auto i = 0; i < 3; ++i) {
    temp1[i] = residual(i);
    temp2[i] = residual(i + 3);
  }
  temp1[0] = -temp1[0];
  temp1[1] = -temp1[1];
  temp1[2] = -temp1[2];
  vReacs[vert1] = temp1;
  vReacs[vert2] = temp2;
}

auto LSFEMBernoulli::getNumberOfIntergrationPoints(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem)
    -> indexType {
  return 2;
}

void LSFEMBernoulli::toParaviewAdaper(PointerCollection &pointers,
                                      FiniteElement::GenericFiniteElement *elem,
                                      vtkPlotInterface &paraviewAdapter,
                                      ParaviewSwitch control)
{
  pointers.getSPDLogger().warn("Plot functionality for element formulation LSFEMBernoulli not implemented!");
}

} // End Namespace

