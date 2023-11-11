// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include <elementFormulations/EL205_HDivTest.h>

#include "plot/vtkplotClassBase.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"


#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include "control/ParameterList.h"
#include <pointercollection/pointercollection.h>



#include "finiteElements/Face.h"

#include <materials/GenericMaterialFormulation.h>

#include <geometry/VertexData.h>

#include <elementFormulations/GenericElementFormulation.h>
#include <solver/GenericSolutionState.h>

#include <MatrixTypes.h>

#include <Eigen/Eigenvalues>
#include <vector>
#include <vtkCellType.h>

namespace HierAMuS::Elementformulations {

EL205_HDivTest::EL205_HDivTest(PointerCollection *ptrCol)
    : GenericElementFormulationInterface(ptrCol) {}

EL205_HDivTest::~EL205_HDivTest() = default;

void EL205_HDivTest::readData(PointerCollection &pointers,
                              ParameterList &list) {

  this->plainstrain = list.getIndexVal("plainstrain");
  this->disporder = list.getIndexVal("disporder");
  this->stressorder = list.getIndexVal("stressorder");
  this->mode = list.getIndexVal("mode");
  this->meshiddisp = list.getIndexVal("meshiddisp");
  this->meshidstress = list.getIndexVal("meshidstress");
  this->E = list.getPrecVal("E");
  this->nu = list.getPrecVal("nu");
}

void EL205_HDivTest::setDegreesOfFreedom(PointerCollection &pointers,
                                         FiniteElement::Face &elem) {

  switch (this->mode) {
  case 1: // Flussproblem
    elem.setHDivShapes(pointers, this->meshidstress, this->stressorder,
                        NodeTypes::displacement);
    break;
  case 2: // hellinger-reissner
    elem.setH1Shapes(pointers, this->meshiddisp, this->disporder);
    elem.setHDivShapes(pointers, this->meshidstress, this->stressorder,
                        NodeTypes::displacement);
    break;
  default:
    throw std::runtime_error("Element formulation 205 called with wrong mode!");
  }
}

void EL205_HDivTest::AdditionalOperations(
  PointerCollection& pointers, FiniteElement::Face &elem) {
  switch (this->mode) {
  case 1: // Flussproblem
    elem.setAllNodeBoundaryConditionMeshId(pointers, this->meshidstress, 2);
    elem.setAllNodeBoundaryConditionMeshId(pointers, this->meshidstress, 1);
    break;
  case 2: // hellinger-reissner
    elem.setAllNodeBoundaryConditionMeshId(pointers, this->meshiddisp, 2);
    elem.setAllNodeBoundaryConditionMeshId(pointers, this->meshidstress, 2);
    break;
  default:
    throw std::runtime_error("Element formulation 205 called with wrong mode!");
  }
}

auto EL205_HDivTest::getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
-> std::vector<DegreeOfFreedom *>  {
  std::vector<DegreeOfFreedom *> Dofs;
  switch (this->mode) {
  case 1: // Flussproblem
    elem->getHDivDofs(pointers, Dofs, this->meshidstress,
                             this->stressorder);
    break;
  case 2: // hellinger-reissner
  {
    std::vector<DegreeOfFreedom *> Dofsstress;

    Dofs.clear();
    elem->getH1Dofs(pointers, Dofs, this->meshiddisp, this->disporder);
    elem->getHDivDofs(pointers, Dofsstress, this->meshidstress,
                      this->stressorder);
    Dofs.insert(Dofs.end(), Dofsstress.begin(), Dofsstress.end());

  }
  break;
  default:
    throw std::runtime_error("Element formulation 205 called with wrong mode!");
  }
  return Dofs;
}

void EL205_HDivTest::setTangentResidual(
  PointerCollection& pointers,
  FiniteElement::Face &elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {
  switch (this->mode) {
  case 1: // Flussproblem
    this->setTangentResidualFluss(pointers, elem, stiffness, residual, Dofs);
    break;
  case 2: // hellinger-reissner
    this->setTangentResidualHellingerReissner(pointers, elem, stiffness,
                                              residual, Dofs);
    break;
  default:
    throw std::runtime_error("Element formulation 205 called with wrong mode!");
  }
}

auto EL205_HDivTest::getNumberOfIntergrationPoints(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem)
    -> indexType  {
  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder((this->stressorder + 1) * 2);
  return GP.getTotalGP();
}

void EL205_HDivTest::setTangentResidualFluss(
  PointerCollection& pointers,
  FiniteElement::Face &elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  Types::VectorX<prec> shapeDeriv;
  Types::Matrix2X<prec> shapeT;
  Types::Matrix2X<prec> coorDeriv;
  Types::VectorXT<prec> coor;
  Types::VectorX<prec> Bmat;
  Types::VectorX<prec> solution;
  std::vector<GenericNodes *> dispNodes;
  Materials::MaterialTransferData materialData;

  Dofs.clear();
  elem.getHDivDofs(pointers, Dofs, this->meshidstress, this->stressorder);

  auto numDofs = static_cast<indexType>(Dofs.size());

  stiffness.resize(numDofs, numDofs);
  stiffness.setZero();

  residual.resize(numDofs);
  residual.setZero();

  elem.getSolution(pointers, Dofs, solution);

  prec da, da2;
  da = 0, da2 = 0;

  Bmat.resize(numDofs);
  Bmat.setZero();

  auto GP = elem.getIntegrationPoints(pointers);
  GP.setOrder((this->stressorder + 1) * 2);

  for (auto &ip : GP) {
    auto jacobi = elem.getJacobian(pointers, ip);
    auto hdivShapes = elem.getHDivShapes(pointers, stressorder, jacobi, ip);

    Bmat = getBMatrixFluss(shapeDeriv, numDofs);
    // materialData.strains = Bmat * solution;

    // elem->getMaterialFormulation()->getMaterialData(materialData);

    da = (jacobi.determinant()) * ip.weight;
    // stiffness += Bmat.transpose() * materialData.materialTangent * Bmat * da;
    stiffness += Bmat * Bmat.transpose() * da;
    da2 += da;
  }
  residual = stiffness * solution;
  
}

Types::VectorX<prec>
EL205_HDivTest::getBMatrixFluss(const Types::VectorX<prec> &shapeDeriv,
                                indexType numDofs) {

  Types::VectorX<prec> Bmat;
  Bmat.resize(numDofs);
  Bmat.setZero();
  for (auto j = 0; j < numDofs / 3; ++j) {
    Bmat(j * 3) = shapeDeriv(j);
  }
  return Bmat;
}

void EL205_HDivTest::setTangentResidualHellingerReissner(
  PointerCollection& pointers,
  FiniteElement::Face &elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  Types::MatrixXX<prec> Bmatstress;
  Types::MatrixXX<prec> Bmatdisp;
  Types::MatrixXX<prec> K;
  Types::VectorX<prec> solution;
  //Materials::MaterialTransferData materialData;
  std::vector<DegreeOfFreedom *> Dofsstress;
  std::vector<DegreeOfFreedom *> Dofsdisp;

  Dofs.clear();
  elem.getH1Dofs(pointers, Dofsdisp, this->meshiddisp, this->disporder);
  elem.getHDivDofs(pointers, Dofsstress, this->meshidstress,
                    this->stressorder);
  Dofs = Dofsdisp;
  Dofs.insert(Dofs.end(), Dofsstress.begin(), Dofsstress.end());

  indexType numDofsstress = Dofsstress.size();
  indexType numDofsdisp = Dofsdisp.size();

  stiffness.resize(numDofsdisp + numDofsstress, numDofsdisp + numDofsstress);
  stiffness.setZero();
  K.resize(numDofsdisp + numDofsstress, numDofsdisp + numDofsstress);
  K.setZero();

  residual.resize(numDofsdisp + numDofsstress);
  residual.setZero();

  elem.getSolution(pointers, Dofs, solution);

  Bmatstress.resize(4, numDofsstress);
  Bmatstress.setZero();
  Bmatdisp.resize(4, numDofsdisp);
  Bmatdisp.setZero();

  prec da = 0;


  auto GP = elem.getIntegrationPoints(pointers);
  GP.setOrder((this->stressorder + 1) * 2);


  Types::MatrixXX<prec> Material = getMaterial(&elem);

  for (auto i: GP) {
    auto jacobi = elem.getJacobian(pointers, i);
    auto H1Shapes = elem.getH1Shapes(pointers, this->disporder, jacobi, i);
    auto HDivShapes =
        elem.getHDivShapes(pointers, this->stressorder, jacobi, i);

    Bmatstress = getBMatrixstress(HDivShapes.shapes, numDofsstress);
    Bmatdisp = getBMatrixepsilon(H1Shapes.shapeDeriv, numDofsdisp);

    da = (jacobi.determinant()) * i.weight;

    // elem->getMaterialFormulation()->getMaterialData(materialData);

    K.topRightCorner(numDofsdisp, numDofsstress) =
        Bmatdisp.transpose() * Bmatstress;
    K.bottomLeftCorner(numDofsstress, numDofsdisp) =
        Bmatstress.transpose() * Bmatdisp;
    K.bottomRightCorner(numDofsstress, numDofsstress) =
        -Bmatstress.transpose() * Material *
        Bmatstress;

    stiffness += K * da;
  }
  residual = stiffness * solution;
}

Types::MatrixXX<prec>
EL205_HDivTest::getBMatrixepsilon(const Types::Matrix2X<prec> &shapeDerivdisp,
                                  indexType numDofsdisp) {

  Types::MatrixXX<prec> Bmatepsilon;
  Bmatepsilon.resize(4, numDofsdisp);
  Bmatepsilon.setZero();
  for (auto j = 0; j < numDofsdisp / 3; ++j) {
    Bmatepsilon(0, j * 3) = shapeDerivdisp(0, j);
    Bmatepsilon(1, j * 3) = 0.5 * shapeDerivdisp(1, j);
    Bmatepsilon(1, j * 3 + 1) = 0.5 * shapeDerivdisp(0, j);
    Bmatepsilon(2, j * 3) = 0.5 * shapeDerivdisp(1, j);
    Bmatepsilon(2, j * 3 + 1) = 0.5 * shapeDerivdisp(0, j);
    Bmatepsilon(3, j * 3 + 1) = shapeDerivdisp(1, j);
  }
  return Bmatepsilon;
}

Types::MatrixXX<prec>
EL205_HDivTest::getBMatrixstress(const Types::Matrix2X<prec> &shapestress,
                                 indexType numDofsstress) {

  Types::MatrixXX<prec> Bmatstress;
  Bmatstress.resize(4, numDofsstress);
  Bmatstress.setZero();
  for (auto j = 0; j < numDofsstress / 3; ++j) {
    Bmatstress(0, j * 3) = shapestress(0, j);
    Bmatstress(1, j * 3) = shapestress(1, j);
    Bmatstress(2, j * 3 + 1) = shapestress(0, j);
    Bmatstress(3, j * 3 + 1) = shapestress(1, j);
  }
  return Bmatstress;
}

Types::MatrixXX<prec>
EL205_HDivTest::getMaterial(FiniteElement::GenericFiniteElement *elem) {

  Types::MatrixXX<prec> Material(4, 4);

  switch (this->plainstrain)
  {
  case 0:
      Material.setZero();
      Material(0, 0) = prec(1.0) / this->E;
      Material(0, 3) = -prec(1.0) / this->E * this->nu;
      Material(1, 1) = (prec(1.0) + this->nu) / (prec(1.0) * this->E);
      Material(2, 2) = (prec(1.0) + this->nu) / (prec(1.0) * this->E);
      Material(3, 0) = -prec(1.0) / this->E * this->nu;
      Material(3, 3) = prec(1.0) / this->E;
    break;
  case 1:
      Material.setZero();
      Material(0, 0) = (prec(1.0) - this->nu * this->nu) / this->E;
      Material(0, 3) = -this->nu * (prec(1.0) + this->nu) / this->E;
      Material(1, 1) = (prec(1.0) + this->nu) / this->E;
      Material(2, 2) = (prec(1.0) + this->nu) / this->E;
      Material(3, 0) = -this->nu * (prec(1.0) + this->nu) / this->E;
      Material(3, 3) = (prec(1.0) - this->nu * this->nu) / this->E;
    break;
  default:
    break;
  }

  return Material;
}

void EL205_HDivTest::toParaviewAdaper(PointerCollection &pointers,
                                      FiniteElement::Face &elem,
                                      vtkPlotInterface &paraviewAdapter,
                                      ParaviewSwitch control) {

  switch (control) {
  case ParaviewSwitch::Mesh: {
    int matNum = static_cast<int>(elem.getMaterial()->getNumber());
    elem.geometryToParaview(pointers, paraviewAdapter, 0, matNum);

    break;
  }
  case ParaviewSwitch::Solution: {
    int matNum = static_cast<int>(elem.getMaterial()->getNumber());

    elem.H1SolutionToParaview(pointers, paraviewAdapter, 0, matNum,
                               this->meshiddisp, this->disporder,
                               paraviewNames::DisplacementName());

    break;
  }
  case ParaviewSwitch::Weights: {
    indexType matNum = elem.getMaterial()->getNumber();

    auto GP = elem.getIntegrationPoints(pointers);
    GP.setOrder((this->stressorder + 1) * 2);

    for (auto i : GP) {

      auto jacobi = elem.getJacobian(pointers, i);

      auto H1Shapes = elem.getH1Shapes(pointers, 1, jacobi, i);

      prec dA = jacobi.determinant() * i.weight;


      for (auto i = 0; i < 4; ++i) {
        auto &Vert = elem.getVertex(pointers, i);
        std::vector<prec> val;
        val.clear();
        val.push_back(H1Shapes.shapes(i) * dA);
        paraviewAdapter.SumPointDataWeighted(0, matNum, val, Vert.getId(), 1,
                                             paraviewNames::weightName());
      }
    }

  } break;
  case ParaviewSwitch::ProjectedValues: {
    indexType matNum = elem.getMaterial()->getNumber();
    auto GP = elem.getIntegrationPoints(pointers);
    switch (this->mode) {
    case 1:
      GP.setOrder((this->stressorder + 1) * 2);
      break;
    case 2:
      GP.setOrder((this->stressorder + 1) * 2);
      break;
    }
    Types::VectorX<prec> solution, solutionStress;
    Materials::MaterialTransferData materialData;

    std::vector<DegreeOfFreedom *> Dofs, DofsStress;

    elem.getH1Dofs(pointers, Dofs, this->meshiddisp, this->disporder);
    elem.getHDivDofs(pointers, DofsStress, this->meshidstress,
                      this->stressorder);
    elem.getSolution(pointers, Dofs, solution);
    elem.getSolution(pointers, DofsStress, solutionStress);

    indexType numDofsstress = DofsStress.size();
    indexType numDofsdisp = Dofs.size();

    Types::MatrixXX<prec> Bmatdisp;
    Bmatdisp.resize(4, Dofs.size());
    Bmatdisp.setZero();
    Types::MatrixXX<prec> Bmatstress;
    Bmatstress.resize(4, DofsStress.size());
    Bmatstress.setZero();

    Types::VectorX<prec> eps, sig, solsig, soleps;
    eps.resize(9);
    sig.resize(9);
    eps.setZero();
    sig.setZero();
    soleps.resize(4);
    solsig.resize(4);
    soleps.setZero();
    solsig.setZero();
    prec da;

    for (auto i:GP) {
      auto jacobi = elem.getJacobian(pointers, i);

      auto H1shapes = elem.getH1Shapes(pointers, this->disporder, jacobi, i);
      auto HDivshapes =
          elem.getHDivShapes(pointers, this->stressorder, jacobi, i);

      Bmatstress = getBMatrixstress(HDivshapes.shapes, numDofsstress);
      Bmatdisp = getBMatrixepsilon(H1shapes.shapeDeriv, numDofsdisp);

      solsig = Bmatstress * solutionStress;
      soleps = Bmatdisp * solution;

      da = (jacobi.determinant()) * i.weight;

      for (auto i = 0; i < 4; ++i) {

        auto &Vert = elem.getVertex(pointers, i);
        std::vector<prec> epsval, sigval;
        epsval.clear();
        sigval.clear();
        prec dAq = da * H1shapes.shapes(i);
        for (auto j = 0; j < 9; ++j) {
          epsval.push_back(prec(0));
          sigval.push_back(prec(0));
        }
        sigval[0] = solsig(0) * dAq;
        sigval[1] = solsig(3) * dAq;
        sigval[3] = solsig(1) * dAq;
        sigval[4] = solsig(2) * dAq;
        sigval[8] = sqrt(solsig(0) * solsig(0) + solsig(3) * solsig(3) - solsig(0) * solsig(3) + 3 * solsig(1) * solsig(2)) * dAq;
        epsval[0] = soleps(0) * dAq;
        epsval[1] = soleps(3) * dAq;
        epsval[3] = soleps(1) * dAq;
        epsval[4] = soleps(2) * dAq;
        paraviewAdapter.SumPointDataWeighted(0, matNum, epsval, Vert.getId(), 9,
                                             paraviewNames::strainName9());
        paraviewAdapter.SumPointDataWeighted(0, matNum, sigval, Vert.getId(), 9,
                                             paraviewNames::stressName9());
      }
    }
  }
  default:
  {}
  }
}

auto  EL205_HDivTest::computeNorm(
        PointerCollection& pointers,
        FiniteElement::Face* elem,
        FiniteElement::NormTypes type) -> prec {
    /*switch (type) {
    case FiniteElement::NormTypes::L2Stresses:*/
  auto GP = elem->getIntegrationPoints(pointers);
    switch (this->mode) {
    case 1:
        GP.setOrder((this->stressorder + 1) * 2);
        break;
    case 2:
        GP.setOrder((this->stressorder + 1) * 2);
        break;
    }
    Types::Matrix2X<prec> shapestress;
    Types::VectorX<prec> shapeDerivstress;
    Types::VectorX<prec> solution, solutionStress;
    std::vector<DegreeOfFreedom*> Dofs, DofsStress;

    prec Norm = 0;
    elem->getHDivDofs(pointers, DofsStress, this->meshidstress,
        this->stressorder);
    elem->getSolution(pointers, DofsStress, solutionStress);

    indexType numDofsstress = DofsStress.size();

    Types::MatrixXX<prec> Bmatstress;
    Bmatstress.resize(4, DofsStress.size());
    Bmatstress.setZero();


    elem->getH1Dofs(pointers, Dofs, this->meshiddisp, this->disporder);
    elem->getSolution(pointers, Dofs, solution);
    //indexType numDofsdisp = Dofs.size();

    Types::MatrixXX<prec> Bmatdisp;
    Bmatdisp.resize(4, Dofs.size());
    Bmatdisp.setZero();

    Types::VectorX<prec> solsig;
    solsig.resize(4);
    solsig.setZero();
    prec da;

    for (auto i:GP) {
      auto jacobi = elem->getJacobian(pointers, i);
      auto HDivShapes =
          elem->getHDivShapes(pointers, this->stressorder, jacobi, i);
        Bmatstress = getBMatrixstress(HDivShapes.shapes, numDofsstress);
        da = (jacobi.determinant()) * i.weight;
        solsig = Bmatstress * solutionStress;
        Norm = Norm + (solsig[1] - solsig[2]) * (solsig[1] - solsig[2]) * da;
    }

    return Norm;
    /*break;
    } */

}
} // namespace HierAMuS
