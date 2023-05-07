// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "forwarddeclaration.h"
#include "plot/vtkplotClassBase.h"



#include <equations/DegreeOfFreedom.h>


#include <elementFormulations/EL102_Timoshenko2D.h>
#include <elementFormulations/GenericElementFormulation.h>

#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

#include <finiteElements/GenericFiniteElement.h>

#include <geometry/Base.h>
#include <geometry/Vertex.h>

#include <equations/EquationHandler.h>
#include <equations/GenericNodes.h>
#include <equations/NodeSet.h>
#include <equations/Nodetypes.h>


#include <pointercollection/pointercollection.h>
#include <solver/GenericSolutionState.h>

#include <Eigen/Dense>

#include <loads/PropfunctionHandler.h>

#include <stdexcept>
#include <types/MatrixTypes.h>
#include <vector>


#include <vtkCellType.h>


namespace HierAMuS::Elementformulations {

EL102_Timoshenko2D::EL102_Timoshenko2D(PointerCollection *ptrCol)
    : GenericElementFormulation(ptrCol) {}

EL102_Timoshenko2D::~EL102_Timoshenko2D() = default;

void EL102_Timoshenko2D::readData(PointerCollection &pointers,
                                  ParameterList &list) {

  this->EA = list.getPrecVal("ea");
  this->EI = list.getPrecVal("ei");
  this->GA = list.getPrecVal("ga");
  this->rhoA = list.getPrecVal("rhoa");

  this->qx = list.getPrecVal("qx");
  this->qy = list.getPrecVal("qy");
  this->mz = list.getPrecVal("mz");


  this->propnum = list.getIndexVal("propnum");
  this->localLoad = list.getIndexVal("local");


  this->meshIdDisp = list.getIndexVal("meshiddisp");
  this->meshIdRot = list.getIndexVal("meshidrot");

  this->disporder = list.getIndexVal("disporder");
  if(this->disporder==0)
    this->disporder=1;
  this->rotorder = list.getIndexVal("rotorder");
  if(this->rotorder==0)
    this->rotorder=1;
  this->mode = list.getIndexVal("mode");
  if(this->mode==0)
    this->mode=1;


  this->messageUnprocessed(pointers, list, "EL102_Timoshenko2D");
}

void EL102_Timoshenko2D::setDegreesOfFreedom(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) {

  if (this->mode == 1) {
    elem->setH1Shapes(pointers, this->meshIdDisp, this->disporder);
    elem->setH1Shapes(pointers, this->meshIdRot, this->rotorder);
  } else {
    elem->setH1Shapes(pointers, this->meshIdDisp, this->disporder);
  }
}

void EL102_Timoshenko2D::AdditionalOperations(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) {
  if (this->mode == 1) {
    elem->setAllNodeBoundaryConditionMeshId(pointers, this->meshIdDisp, 2);
    elem->setAllNodeBoundaryConditionMeshId(pointers, this->meshIdRot, 1);
    elem->setAllNodeBoundaryConditionMeshId(pointers, this->meshIdRot, 2);
  }
}

auto EL102_Timoshenko2D::getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
-> std::vector<DegreeOfFreedom *>  {

  std::vector<DegreeOfFreedom *> Dofs, dispDofs, rotDofs;
  elem->getH1Dofs(pointers, dispDofs, this->meshIdDisp, this->disporder);
  elem->getH1Dofs(pointers, rotDofs, this->meshIdRot, this->rotorder);
  Dofs.clear();
  Dofs.insert(Dofs.end(), dispDofs.begin(), dispDofs.end());
  Dofs.insert(Dofs.end(), rotDofs.begin(), rotDofs.end());
  return Dofs;
}

void EL102_Timoshenko2D::setTangentResidual(
  PointerCollection& pointers,
  FiniteElement::GenericFiniteElement *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  switch (this->mode) {
  case 1: {
    this->setTangentResidualLinear(pointers, elem, stiffness, residual, Dofs);
  } break;
  case 2: {
    this->setTangentResidualNonlinear(pointers, elem, stiffness, residual,
                                      Dofs);
  } break;
  default: {
    throw std::runtime_error("\nSelected mode for EL102 not implemented!");
  }
  }
}

void EL102_Timoshenko2D::setTangentResidualLinear(
  PointerCollection& pointers,
  FiniteElement::GenericFiniteElement *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {
  
  indexType intorder;
  this->disporder > this->rotorder ? intorder = this->disporder
                                   : intorder = this->rotorder;
  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(intorder*2);
  

  
  std::vector<DegreeOfFreedom *> dispDofs, rotDofs;
  elem->getH1Dofs(pointers, dispDofs, this->meshIdDisp, this->disporder);
  elem->getH1Dofs(pointers, rotDofs, this->meshIdRot, this->rotorder);
  Dofs.clear();
  Dofs.insert(Dofs.end(), dispDofs.begin(), dispDofs.end());
  Dofs.insert(Dofs.end(), rotDofs.begin(), rotDofs.end());

  Types::VectorX<prec> disp;
  elem->getSolution(pointers, Dofs, disp);
    

  stiffness.resize(Dofs.size(), Dofs.size());
  residual.resize(Dofs.size());
  stiffness.setZero();
  residual.setZero();

  indexType numDispNodes = dispDofs.size() / 3;
  indexType numRotNodes = rotDofs.size() / 3;

  Types::Matrix3X<prec> Bmat;

  Bmat.resize(3, Dofs.size());
  Bmat.setZero();
  

  Types::Matrix33<prec> cmat = this->getMaterialMatrix();


  for (auto i : GP) {
    Types::Vector3<prec> dirVec;
    Types::Vector3<prec> dirVec2;
    Types::MatrixXX<prec> jacobi;
    dirVec = elem->getA1Vector(pointers, i);
    dirVec2(0) = -dirVec(1);
    dirVec2(1) = dirVec(0);
    dirVec2(2) = prec(0);
    jacobi = elem->getJacobian(pointers, i);
    auto dispShapes = elem->getH1Shapes(pointers, this->disporder, jacobi, i);
    auto rotShapes = elem->getH1Shapes(pointers, this->rotorder, jacobi, i);
    

    for (auto i = 0; i < numDispNodes; ++i) {
      Bmat.block(0, 3 * i, 1, 3) = dirVec.transpose() * dispShapes.shapeDeriv(i);
      Bmat.block(1, 3 * i, 1, 3) =
          dirVec2.transpose() * dispShapes.shapeDeriv(i);
    }

    indexType offset = 3 * numDispNodes;
    for (auto i = 0; i < numRotNodes; ++i) {
      Bmat(2, 3 * i + offset) = -rotShapes.shapeDeriv(i);
      Bmat(1, 3 * i + offset) = -rotShapes.shapes(i);
    }

    stiffness += Bmat.transpose() * cmat * Bmat * jacobi.determinant() * i.weight;
  }

  residual = stiffness * disp;
}

void EL102_Timoshenko2D::setTangentResidualNonlinear(
  PointerCollection& pointers,
  FiniteElement::GenericFiniteElement *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  auto &vert1 = elem->getVertex(pointers, 0);
  auto &vert2 = elem->getVertex(pointers, 1);

  Types::Vector3<prec> coor1, coor2;
  coor1 = vert1.getCoordinates();
  coor2 = vert2.getCoordinates();
  prec length, detj;
  elem->getJacobian(pointers, detj, 0);
  length = (prec)2 * detj;

  Dofs.clear();
  elem->getH1Dofs(pointers, Dofs, this->meshIdDisp, this->disporder);
  auto sol = pointers.getSolutionState();
  Eigen::Matrix<prec, Eigen::Dynamic, 1> disp;
  disp = sol->getSolution(Dofs);

  // std::cout << disp.transpose() << std::endl;

  auto matsizes = static_cast<indexType>(Dofs.size());
  indexType numNodes = matsizes / 3;

  std::vector<prec> gp, weight;
  elem->getGaussPoints(this->disporder + 1, weight, gp);

  Eigen::Matrix<prec, 2, 1> svec;
  svec(0) = coor2[0] - coor1[0];
  svec(1) = coor2[1] - coor1[1];

  svec = svec / length;

  prec css, sss;

  css = svec(0);
  sss = svec(1);

  Types::MatrixXX<prec> T;
  T.resize(matsizes, matsizes);
  T.setZero();
  for (auto i = 0; i < numNodes; ++i) {
    T(i * 3, i * 3) = css;
    T(i * 3, i * 3 + 1) = sss;
    T(i * 3 + 1, i * 3) = -sss;
    T(i * 3 + 1, i * 3 + 1) = css;
    T(i * 3 + 2, i * 3 + 2) = (prec)1.0;
  }

  Eigen::Matrix<prec, Eigen::Dynamic, 1> localDisp;
  localDisp = T * disp;

  Eigen::Matrix<prec, 3, Eigen::Dynamic> Bmat;
  Eigen::Matrix<prec, 3, 1> strains, stress;
  Eigen::Matrix<prec, 3, 3> mat;

  Types::MatrixXX<prec> GN, GQ, GM;
  mat.setZero();
  mat(0, 0) = this->EA;
  mat(1, 1) = this->GA;
  mat(2, 2) = this->EI;
  Bmat.resize(3, matsizes);
  Bmat.setZero();
  stiffness.resize(matsizes, matsizes);
  stiffness.setZero();
  residual.resize(matsizes);
  residual.setZero();
  GN.resize(matsizes, matsizes);
  GQ.resize(matsizes, matsizes);
  GM.resize(matsizes, matsizes);
  for (auto i = 0; i < gp.size(); ++i) {
    Types::VectorX<prec> shape, shapeDeriv;
    elem->getH1Shapes(pointers, this->disporder, detj, shape, shapeDeriv,
                      gp[i]);

    prec ucx = 0, wcx = 0, bcx = 0, beta = 0;

    for (auto i = 0; i < numNodes; ++i) {
      ucx += shapeDeriv(i) * localDisp(3 * i);
      wcx += shapeDeriv(i) * localDisp(3 * i + 1);
      bcx += shapeDeriv(i) * localDisp(3 * i + 2);
      beta += shape(i) * localDisp(3 * i + 2);
    }
    // prec ucx = shapeDeriv(0) * localDisp(0) + shapeDeriv(1) * localDisp(3);
    // prec wcx = shapeDeriv(0) * localDisp(1) + shapeDeriv(1) * localDisp(4);
    // prec bcx = shapeDeriv(0) * localDisp(2) + shapeDeriv(1) * localDisp(5);
    // prec beta = shape(0) * localDisp(2) + shape(1) * localDisp(5);

    strains(0) = ucx + (prec)0.5 * ucx * ucx + (prec)0.5 * wcx * wcx; // epsilon
    strains(1) = ((prec)1 + ucx) * sin(beta) + wcx * cos(beta);       // gamma
    strains(2) = (((prec)1 + ucx) * cos(beta) - wcx * sin(beta)) * bcx; // kappa

    stress = mat * strains;
    prec alpha1 = (-((prec)1 + ucx) * sin(beta) - wcx * cos(beta));
    prec alpha2 = (((prec)1 + ucx) * cos(beta) - wcx * sin(beta));
    for (auto j = 0; j < numNodes; ++j) {
      Bmat(0, j * 3) = ((prec)1.0 + ucx) * shapeDeriv(j);
      Bmat(0, j * 3 + 1) = wcx * shapeDeriv(j);

      Bmat(1, j * 3) = sin(beta) * shapeDeriv(j);
      Bmat(1, j * 3 + 1) = cos(beta) * shapeDeriv(j);
      Bmat(1, j * 3 + 2) = alpha2 * shape(j);

      Bmat(2, j * 3) = cos(beta) * bcx * shapeDeriv(j);
      Bmat(2, j * 3 + 1) = -bcx * sin(beta) * shapeDeriv(j);
      Bmat(2, j * 3 + 2) = alpha1 * bcx * shape(i) + alpha2 * shapeDeriv(j);
    }
    GN.setZero();
    GQ.setZero();
    GM.setZero();
    for (auto j = 0; j < numNodes; ++j) {
      for (auto k = 0; k < numNodes; ++k) {
        GN(j * 3, k * 3) = shapeDeriv(j) * shapeDeriv(k);
        GN(j * 3 + 1, k * 3 + 1) = shapeDeriv(j) * shapeDeriv(k);

        GQ(j * 3, k * 3 + 2) = shapeDeriv(j) * cos(beta) * shape(k);
        GQ(j * 3 + 1, k * 3 + 2) = -shapeDeriv(j) * sin(beta) * shape(k);
        GQ(j * 3 + 2, k * 3) = shape(j) * cos(beta) * shapeDeriv(k);
        GQ(j * 3 + 2, k * 3 + 1) = -shape(j) * sin(beta) * shapeDeriv(k);
        GQ(j * 3 + 2, k * 3 + 2) = shape(j) * alpha1 * shape(k);

        prec a = -shape(j) * sin(beta) * bcx * shapeDeriv(k) +
                 shapeDeriv(j) * cos(beta) * shapeDeriv(k);
        prec b = -shape(j) * cos(beta) * bcx * shapeDeriv(k) -
                 shapeDeriv(j) * sin(beta) * shapeDeriv(k);
        prec c = -shapeDeriv(j) * sin(beta) * bcx * shape(k) +
                 shapeDeriv(j) * cos(beta) * shapeDeriv(k);
        prec d = -shapeDeriv(j) * cos(beta) * bcx * shape(k) -
                 shapeDeriv(j) * sin(beta) * shapeDeriv(k);
        prec e = -shape(j) * alpha2 * bcx * shape(k) +
                 shape(j) * alpha1 * shapeDeriv(k) +
                 shapeDeriv(j) * alpha1 * shape(k);
        GM(j * 3, k * 3 + 2) = c;
        GM(j * 3 + 1, k * 3 + 2) = d;
        GM(j * 3 + 2, k * 3 + 2) = e;
        GM(j * 3 + 2, k * 3) = a;
        GM(j * 3 + 2, k * 3 + 1) = b;
      }
    }

    GN *= stress(0);
    GQ *= stress(1);
    GM *= stress(2);

    residual += Bmat.transpose() * stress * detj * weight[i];
    stiffness +=
        (Bmat.transpose() * mat * Bmat + GN + GQ + GM) * detj * weight[i];
  }
  stiffness = T.transpose() * stiffness * T;
  residual = T.transpose() * residual;

  prec propVal;
  auto tempProp = pointers.getPropLoads();

  propVal = tempProp->getPropValue(this->propnum);

  residual(0) -= this->qx * propVal / (prec)2 * length;
  residual(3) -= this->qx * propVal / (prec)2 * length;
  residual(1) -= this->qy * propVal / (prec)2 * length;
  residual(4) -= this->qy * propVal / (prec)2 * length;
}

void EL102_Timoshenko2D::setMass(
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
  // set2->getNodes(&this->ptrCol,temp);
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

  Types::Vector2<prec> svec;
  svec(0) = coor2(0) - coor1(0);
  svec(1) = coor2(1) - coor1(1);

  svec = svec / length;

  prec css, sss;

  css = svec(0);
  sss = svec(1);

  auto matsizes = static_cast<indexType>(Dofs.size());
  indexType numNodes = matsizes / 3;

  Types::MatrixXX<prec> T;
  T.resize(matsizes, matsizes);
  T.setZero();
  for (auto i = 0; i < numNodes; ++i) {
    T(i * 3, i * 3) = css;
    T(i * 3, i * 3 + 1) = sss;
    T(i * 3 + 1, i * 3) = -sss;
    T(i * 3 + 1, i * 3 + 1) = css;
    T(i * 3 + 2, i * 3 + 2) = (prec)1.0;
  }

  stiffness.resize(6, 6);
  for (auto i = 0; i < 2; ++i) {
    stiffness(i, i) = this->rhoA * length / (prec)2;
    stiffness(i + 3, i + 3) = this->rhoA * length / (prec)2;
  }
  stiffness(5, 5) = this->rhoA * length / (prec)2 * (prec)1e-6;
  stiffness(2, 2) = this->rhoA * length / (prec)2 * (prec)1e-6;

  stiffness = T.transpose() * stiffness * T;
}

void EL102_Timoshenko2D::getElementsLocalNodalReactions(
    PointerCollection &ptrCol, FiniteElement::GenericFiniteElement *elem,
    std::map<indexType, std::vector<prec>> &vReacs) {
  // indexType vert1 = elem->getVertexId(ptrCol, 0);
  // indexType vert2 = elem->getVertexId(ptrCol, 1);

  // Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> stiffness;
  // Eigen::Matrix<prec, Eigen::Dynamic, 1> residual;
  // std::vector<DegreeOfFreedom *> Dofs;

  // this->setTangentResidual(elem, stiffness, residual, Dofs);

  // GenericGeometryElement *gvert1, *gvert2;
  // gvert1 = elem->getVertex(*this->ptrCol, 0);
  // gvert2 = elem->getVertex(*this->ptrCol, 1);

  // std::vector<prec> coor1, coor2;
  // coor1 = gvert1->getCoordinates();
  // coor2 = gvert2->getCoordinates();
  // prec length = (prec)0, tempDiff;
  // for (auto i = 0; i < 3; ++i) {
  //   tempDiff = (coor1[i] - coor2[i]);
  //   length += tempDiff * tempDiff;
  // }
  // length = sqrt(length);

  // Eigen::Matrix<prec, 2, 1> svec;
  // svec(0) = coor2[0] - coor1[0];
  // svec(1) = coor2[1] - coor1[1];

  // svec = svec / length;

  // prec css, sss;

  // css = svec(0);
  // sss = svec(1);

  // Eigen::Matrix<prec, 6, 6> T;
  // T.setZero();
  // T(0, 0) = css;
  // T(0, 1) = sss;
  // T(1, 0) = -sss;
  // T(1, 1) = css;
  // T(2, 2) = (prec)1.0;

  // T(3, 3) = css;
  // T(3, 4) = sss;
  // T(4, 3) = -sss;
  // T(4, 4) = css;
  // T(5, 5) = (prec)1.0;

  // residual = T * residual;

  // std::vector<prec> temp1(6), temp2(6);
  // for (auto i = 0; i < 3; ++i) {
  //   temp1[i] = residual(i);
  //   temp2[i] = residual(i + 3);
  // }
  // temp1[0] = -temp1[0];
  // temp1[1] = -temp1[1];
  // temp1[2] = -temp1[2];
  // vReacs[vert1] = temp1;
  // vReacs[vert2] = temp2;
}

void EL102_Timoshenko2D::toParaviewAdaper(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem,
    vtkPlotInterface &paraviewAdapter, ParaviewSwitch control) {

  switch (control) {
  case HierAMuS::ParaviewSwitch::Mesh: {
    int matNum = static_cast<int>(elem->getMaterial()->getNumber());

    Types::Vector3<prec> coors;
    std::vector<indexType> cell;
    cell.resize(2);
    for (auto i = 0; i < 2; ++i) {
      auto &Vert = elem->getVertex(pointers, i);
      coors = Vert.getCoordinates();
      indexType id = Vert.getId();
      paraviewAdapter.addPoint(0, matNum, id, coors(0), coors(1), coors(2));
      cell[i] = id;
    }
    int vtkNumber = VTK_LINE;
    indexType elid = elem->getId();
    paraviewAdapter.addCell(0, matNum, elid, 1, cell, 2, vtkNumber);

    break;
  }
  case HierAMuS::ParaviewSwitch::Solution: {
    int matNum = static_cast<int>(elem->getMaterial()->getNumber());

    std::vector<DegreeOfFreedom *> Dofs;
    elem->getH1Dofs(pointers, Dofs, this->meshIdDisp, this->disporder);
    Types::VectorX<prec> dispSolution, rotSolution;

    elem->getSolution(pointers, Dofs, dispSolution);

    Dofs.clear();
    elem->getH1Dofs(pointers, Dofs, this->meshIdRot, this->rotorder);
    elem->getSolution(pointers, Dofs, rotSolution);

    indexType cc = 0;
    std::vector<prec> sol(3), solrot(3);
    for (auto i = 0; i < 2; ++i) {
      auto &Vert = elem->getVertex(pointers, i);
      indexType Vid = Vert.getId();
      for (auto j = 0; j < 3; ++j) {
        sol[j] = dispSolution(cc + j);
        solrot[j] = rotSolution(cc + j);
      }
      paraviewAdapter.setPointData(0, matNum, Vid, sol, 3,
                                   paraviewNames::DisplacementName());
      paraviewAdapter.setPointData(0, matNum, Vid, solrot, 3,
                                   paraviewNames::RotationName());

      cc += 3;
    }

    break;
  }
  default:{
    break;
  }
  }

}



auto EL102_Timoshenko2D::getHistoryDataStructure()
    -> const HistoryDataStructure & {
  return m_HistoryDataStructure;
}

auto EL102_Timoshenko2D::getNumberOfIntergrationPoints(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> indexType {
  auto GP = elem->getIntegrationPoints(pointers);
  indexType order = std::max(this->disporder, this->rotorder);
  GP.setOrder(order * 2);
  return GP.getTotalGP();
}

Types::Matrix33<prec> EL102_Timoshenko2D::getMaterialMatrix() const {
  Types::Matrix33<prec> mat;
  mat.setZero();
  mat(0, 0) = this->EA;
  mat(1, 1) = this->GA;
  mat(2, 2) = this->EI;
  return mat;
}

const HistoryDataStructure
    EL102_Timoshenko2D::m_HistoryDataStructure({{4, 2}, {8, 1}}, {});


} // namespace HierAMuS
