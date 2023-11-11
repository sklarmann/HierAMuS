// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause




#include <elementFormulations/EL303_ThermoMechanikSolid3D.h>
#include <pointercollection/pointercollection.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include "control/ParameterList.h"


#include <finiteElements/Volume.h>



#include <materials/GenericMaterialFormulation.h>

#include "shapefunctions/IntegrationsPoints/dataClasses/GaussPoints.h"

#include <elementFormulations/GenericElementFormulation.h>
#include <solver/GenericSolutionState.h>

#include <types/MatrixTypes.h>

#include <Eigen/Eigenvalues>

namespace HierAMuS {
namespace Elementformulations {

EL303_ThermoMechanikSolid3D::EL303_ThermoMechanikSolid3D(PointerCollection *ptrCol)
    : GenericElementFormulationInterface(ptrCol) {}

EL303_ThermoMechanikSolid3D::~EL303_ThermoMechanikSolid3D() {}

void EL303_ThermoMechanikSolid3D::readData(PointerCollection &pointers,
                                           ParameterList &list) {
  this->meshIdDisp = list.getIndexVal("meshiddisp");
  this->meshidTemp = list.getIndexVal("meshidtemperatur");
  this->intOrderDisp = list.getIndexVal("disporder");

  this->mu = list.getPrecVal("mu");
  this->lambda = list.getPrecVal("lambda");
  this->alpha = list.getPrecVal("alpha");
  this->c = list.getPrecVal("c");
  this->rho0 = list.getPrecVal("rho0");
  this->T0 = list.getPrecVal("T0");
  this->kappa = list.getPrecVal("kappa");

  this->mode = list.getIndexVal("mode");
  

  auto &Log = pointers.getSPDLogger();

  Log.info("\n{:-<100}\n"
                "*   Element 303, specified Options\n"
                "    Mesh id for displacement nodes:         {:>12}\n"
                "    Shape function order for displacements: {:>12}\n"
                "    Material parameters:\n"
                "      mu:                                   {:>12.4e}\n"
                "      lambda:                               {:>12.4e}\n"
                "      alpha:                                {:>12.4e}\n"
                "      c:                                    {:>12.4e}\n"
                "      rho0:                                 {:>12.4e}\n"
                "      T0:                                   {:>12.4e}\n"
                "      kappa:                                {:>12.4e}\n"
                "{:-<100}\n",
                "",
                this->meshIdDisp,
                this->intOrderDisp,
                this->mu,
                this->lambda,
                this->alpha,
                this->c,
                this->rho0,
                this->T0,
                this->kappa,
                "");


  this->messageUnprocessed(pointers, list, "EL303_ThermoMechanikSolid3D");
}

void EL303_ThermoMechanikSolid3D::setDegreesOfFreedom(
    PointerCollection &pointers, FiniteElement::Volume &elem) {

  elem.setH1Shapes(pointers, this->meshIdDisp, this->intOrderDisp);
  elem.setH1Shapes(pointers, this->meshidTemp, this->intOrderDisp);
}

void EL303_ThermoMechanikSolid3D::AdditionalOperations(
  PointerCollection& pointers, FiniteElement::Volume &elem) {
  auto vol = elem.getVolume(pointers, 0);
  vol->setAllNodeBoundaryConditionMeshId(this->meshidTemp, 1);
  vol->setAllNodeBoundaryConditionMeshId(this->meshidTemp, 2);

}

auto EL303_ThermoMechanikSolid3D::getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> std::vector<DegreeOfFreedom*>{
  std::vector<DegreeOfFreedom *> ret, temp;

  elem->getH1Dofs(pointers, ret, this->meshIdDisp, this->intOrderDisp);
  elem->getH1Dofs(pointers, temp, this->meshidTemp, this->intOrderDisp);
  ret.insert(ret.end(), temp.begin(), temp.end());
  return ret;
}

void EL303_ThermoMechanikSolid3D::setTangentResidual(
  PointerCollection& pointers, FiniteElement::Volume &elem,
  Types::MatrixXX<prec> &stiffness, Types::VectorX<prec> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  switch (this->mode) {
  case 1:
    this->setTangentResidualLinear(pointers, elem, stiffness, residual, Dofs);
    break;
  case 2:
    this->setTangentResidualNonLinear(pointers, elem, stiffness, residual,
                                      Dofs);
    break;
  }
}

void EL303_ThermoMechanikSolid3D::setTangentResidualLinear(
  PointerCollection& pointers,
  FiniteElement::Volume &elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  Types::VectorX<prec> solution, shapeMat;
  Types::Matrix3X<prec> tempGradient;



  Dofs.clear();
  Dofs = this->getDofs(pointers, &elem);

  indexType numDofs = static_cast<indexType>(Dofs.size());
  indexType numDispDofs = numDofs / 2;

  shapeMat.resize(numDispDofs);
  shapeMat.setZero();
  tempGradient.resize(3, numDispDofs);
  tempGradient.setZero();

  stiffness.resize(numDofs, numDofs);
  stiffness.setZero();

  residual.resize(numDofs);
  residual.setZero();

  elem.getSolution(pointers, Dofs, solution);


  auto GP = elem.getIntegrationPoints(pointers);
  GP.setOrder(this->intOrderDisp * 2);

  auto Hist = elem.getHistoryDataIterator(pointers);

  Materials::MaterialTransferData materialData;
  materialData.historyData = &Hist;
  //indexType cc = 0;
  for (auto i : GP) {
    auto jaco = elem.getJacobian(pointers, i);
    auto shapes = elem.getH1Shapes(pointers, this->intOrderDisp, jaco, i);

    auto Bmat = this->getLinearBMatrix(shapes);
    auto TrMat = this->getLinearBMatrixTrace(shapes);
    for (auto j=0;j<numDispDofs/3;++j)
    {
      shapeMat(3 * j) = shapes.shapes(j);
      tempGradient(0, 3 * j) = shapes.shapeDeriv(0, j);
      tempGradient(1, 3 * j) = shapes.shapeDeriv(1, j);
      tempGradient(2, 3 * j) = shapes.shapeDeriv(2, j);
    }
    //materialData.strains = Bmat * solution;
    //elem->getMaterialFormulation()->getMaterialData(materialData);

    auto dV = i.weight * jaco.determinant();
    // Mechanical part seems to be ok
    stiffness.block(0, 0, numDispDofs, numDispDofs) +=
        Bmat.transpose() * Bmat * dV * this->mu * prec(2);
    stiffness.block(0, 0, numDispDofs, numDispDofs) +=
        TrMat * TrMat.transpose() * dV * this->lambda;

    // Thermal mechanical coupling
    stiffness.block(0, numDispDofs, numDispDofs, numDispDofs) -=
        TrMat * shapeMat.transpose() * dV * this->alpha *
        (prec(3) * this->lambda + prec(2) * this->mu);

    // Thermal part, something seems wrong
    //stiffness.block(numDispDofs, numDispDofs, numDispDofs, numDispDofs) -=
    //    shapeMat * shapeMat.transpose() * this->rho0 / this->T0 * this->c * dV;
    stiffness.block(numDispDofs, numDispDofs, numDispDofs, numDispDofs) -=
        tempGradient.transpose() * tempGradient * this->kappa * dV;

    residual.block(0, 0, numDispDofs, 1) -=
        TrMat * this->T0 * this->alpha* (prec(3)*this->lambda + prec(2)*this->mu) * dV;
    //residual.block(numDispDofs, 0, numDispDofs, 1) +=
    //    shapeMat * this->T0 * this->rho0*this->c/this->T0* dV;

    Hist.next();
  }
  //stiffness.block(numDispDofs, 0, numDispDofs, numDispDofs) =
  //stiffness.block(0, numDispDofs, numDispDofs, numDispDofs).transpose();

  residual += stiffness * solution;

  //std::cout << stiffness.eigenvalues() << std::endl;
  if (elem.getId()==-1) {
    std::cout << stiffness.block(0, 0, numDispDofs, numDispDofs) << "\n"
              << std::endl;
    std::cout << stiffness.block(numDispDofs, 0, numDispDofs, numDispDofs)
              << "\n"
              << std::endl;
    std::cout << stiffness.block(numDispDofs, numDispDofs, numDispDofs, numDispDofs)
              << "\n"
              << std::endl;
  }
}

void EL303_ThermoMechanikSolid3D::setTangentResidualNonLinear(
  PointerCollection& pointers,
  FiniteElement::Volume &elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  Types::VectorX<prec> solution, shapeMat;
  Types::Matrix3X<prec> tempGradient;



  Dofs.clear();
  Dofs = this->getDofs(pointers, &elem);

  indexType numDofs = static_cast<indexType>(Dofs.size());
  indexType numDispDofs = numDofs / 2;

  shapeMat.resize(numDispDofs);
  shapeMat.setZero();
  tempGradient.resize(3, numDispDofs);
  tempGradient.setZero();

  stiffness.resize(numDofs, numDofs);
  stiffness.setZero();

  residual.resize(numDofs);
  residual.setZero();

  elem.getSolution(pointers, Dofs, solution);


  auto GP = elem.getIntegrationPoints(pointers);
  GP.setOrder(this->intOrderDisp * 2);

  auto Hist = elem.getHistoryDataIterator(pointers);

  Materials::MaterialTransferData materialData;
  materialData.historyData = &Hist;
  //indexType cc = 0;
  for (auto i : GP) {
    auto jaco = elem.getJacobian(pointers, i);
    auto shapes = elem.getH1Shapes(pointers, this->intOrderDisp, jaco, i);

    auto Bmat = this->getLinearBMatrix(shapes);
    auto TrMat = this->getLinearBMatrixTrace(shapes);
    for (auto j=0;j<numDispDofs/3;++j)
    {
      shapeMat(3 * j) = shapes.shapes(j);
      tempGradient(0, 3 * j) = shapes.shapeDeriv(0, j);
      tempGradient(1, 3 * j) = shapes.shapeDeriv(1, j);
      tempGradient(2, 3 * j) = shapes.shapeDeriv(2, j);
    }
    for (auto j=0;j<numDispDofs/3;++j)
    {
      shapeMat(3 * j) = shapes.shapes(j);
      tempGradient(0, 3 * j) = shapes.shapeDeriv(0, j);
      tempGradient(1, 3 * j) = shapes.shapeDeriv(1, j);
      tempGradient(2, 3 * j) = shapes.shapeDeriv(2, j);
    }
    //materialData.strains = Bmat * solution;
    //elem->getMaterialFormulation()->getMaterialData(materialData);

    auto dV = i.weight * jaco.determinant();

    // Mechanical part seems to be ok
    stiffness.block(0, 0, numDispDofs, numDispDofs) +=
        Bmat.transpose() * Bmat * dV * this->mu * prec(2);
    stiffness.block(0, 0, numDispDofs, numDispDofs) +=
        TrMat * TrMat.transpose() * dV * this->lambda;

    // Thermal mechanical coupling
    stiffness.block(0, numDispDofs, numDispDofs, numDispDofs) -=
        TrMat * shapeMat.transpose() * dV * this->alpha *
        (prec(3) * this->lambda + prec(2) * this->mu);

    // Thermal part, something seems wrong
    //stiffness.block(numDispDofs, numDispDofs, numDispDofs, numDispDofs) -=
    //    shapeMat * shapeMat.transpose() * this->rho0 / this->T0 * this->c * dV;
    stiffness.block(numDispDofs, numDispDofs, numDispDofs, numDispDofs) -=
        tempGradient.transpose() * tempGradient * this->kappa * dV;

    residual.block(0, 0, numDispDofs, 1) -=
        TrMat * this->T0 * this->alpha* (prec(3)*this->lambda + prec(2)*this->mu) * dV;
    //residual.block(numDispDofs, 0, numDispDofs, 1) +=
    //    shapeMat * this->T0 * this->rho0*this->c/this->T0* dV;

    Hist.next();
  }
  //stiffness.block(numDispDofs, 0, numDispDofs, numDispDofs) =
  //stiffness.block(0, numDispDofs, numDispDofs, numDispDofs).transpose();

  residual += stiffness * solution;

  //std::cout << stiffness.eigenvalues() << std::endl;
  if (elem.getId()==-1) {
    std::cout << stiffness.block(0, 0, numDispDofs, numDispDofs) << "\n"
              << std::endl;
    std::cout << stiffness.block(numDispDofs, 0, numDispDofs, numDispDofs)
              << "\n"
              << std::endl;
    std::cout << stiffness.block(numDispDofs, numDispDofs, numDispDofs, numDispDofs)
              << "\n"
              << std::endl;
  }
}

auto EL303_ThermoMechanikSolid3D::getDeformationGradient(Geometry::H1Shapes &shapes,
                                                 Types::VectorX<prec> &solution)
    -> Types::Matrix33<prec> {
  Types::Matrix33<prec> F;
  F.setIdentity();

  for (auto i = 0; i < shapes.shapes.rows(); ++i) {
    for (auto j = 0; j < 3; ++j) {
      F(0, j) += solution(i * 3) * shapes.shapeDeriv(j, i);
      F(1, j) += solution(i * 3 + 1) * shapes.shapeDeriv(j, i);
      F(2, j) += solution(i * 3 + 2) * shapes.shapeDeriv(j, i);
    }
  }
  return F;
}

auto EL303_ThermoMechanikSolid3D::getNonLinearBMatrix(Geometry::H1Shapes &shapes,
                                              Types::Matrix33<prec> &F)
    -> Types::Matrix6X<prec> {
  indexType nshapes = shapes.shapes.rows();
  Types::Matrix6X<prec> Bmat(6,nshapes*3);
  //Bmat.resize(6, nshapes * 3);
  Bmat.setZero();
  Types::Matrix33<prec> FT = F.transpose();
  auto dx1 = FT.block(0, 0, 1, 3);
  auto dx2 = FT.block(1, 0, 1, 3);
  auto dx3 = FT.block(2, 0, 1, 3);
  for (auto i = 0; i < nshapes; ++i) {
    Bmat.block(0, 3 * i, 1, 3) = shapes.shapeDeriv(0, i) * dx1;
    Bmat.block(1, 3 * i, 1, 3) = shapes.shapeDeriv(1, i) * dx2;
    Bmat.block(2, 3 * i, 1, 3) = shapes.shapeDeriv(2, i) * dx3;

    Bmat.block(3, 3 * i, 1, 3) =
        (shapes.shapeDeriv(0, i) * dx2 + shapes.shapeDeriv(1, i) * dx1);
    Bmat.block(4, 3 * i, 1, 3) =
        (shapes.shapeDeriv(0, i) * dx3 + shapes.shapeDeriv(2, i) * dx1);
    Bmat.block(5, 3 * i, 1, 3) =
        (shapes.shapeDeriv(1, i) * dx3 + shapes.shapeDeriv(2, i) * dx2);
  }

  return Bmat;
}

auto EL303_ThermoMechanikSolid3D::getLinearBMatrix(Geometry::H1Shapes &shapes)
    -> Types::Matrix6X<prec> {
  Types::Matrix6X<prec> Bmat;
  indexType nshapes = shapes.shapes.rows();
  Bmat.resize(6, nshapes * 3);
  Bmat.setZero();
  for (auto l = 0; l < nshapes; ++l) {
    Bmat(0, l * 3) = shapes.shapeDeriv(0, l);
    Bmat(1, l * 3 + 1) = shapes.shapeDeriv(1, l);
    Bmat(2, l * 3 + 2) = shapes.shapeDeriv(2, l);

    Bmat(3, l * 3) = shapes.shapeDeriv(1, l) / sqrt(prec(2));
    Bmat(3, l * 3 + 1) = shapes.shapeDeriv(0, l) / sqrt(prec(2));

    Bmat(4, l * 3) = shapes.shapeDeriv(2, l) / sqrt(prec(2));
    Bmat(4, l * 3 + 2) = shapes.shapeDeriv(0, l) / sqrt(prec(2));

    Bmat(5, l * 3 + 1) = shapes.shapeDeriv(2, l) / sqrt(prec(2));
    Bmat(5, l * 3 + 2) = shapes.shapeDeriv(1, l) / sqrt(prec(2));
  }

  return Bmat;
}

auto EL303_ThermoMechanikSolid3D::getLinearBMatrixTrace(
    Geometry::H1Shapes &shapes) -> Types::VectorX<prec> {

  indexType numshapes = shapes.shapes.size();

  Types::VectorX<prec> ret;
  ret.resize(numshapes * 3);
  ret.setZero();

  for (auto i=0;i<numshapes;++i) {
    ret(3 * i + 0) = shapes.shapeDeriv(0, i);
    ret(3 * i + 1) = shapes.shapeDeriv(1, i);
    ret(3 * i + 2) = shapes.shapeDeriv(2, i);
  }
  return ret;
}

auto EL303_ThermoMechanikSolid3D::getGeometricMatrix(Geometry::H1Shapes &shapes,
                                            Types::VectorX<prec> &stress)
    -> Types::MatrixXX<prec> {
  Types::MatrixXX<prec> G;
  indexType nshapes = shapes.shapes.rows();
  G.resize(nshapes * 3, nshapes * 3);
  G.setZero();

  for (auto i = 0; i < nshapes; ++i) {
    for (auto j = 0; j < nshapes; ++j) {
      prec gik;
      gik = stress(0) * shapes.shapeDeriv(0, i) * shapes.shapeDeriv(0, j);
      gik += stress(1) * shapes.shapeDeriv(1, i) * shapes.shapeDeriv(1, j);
      gik += stress(2) * shapes.shapeDeriv(2, i) * shapes.shapeDeriv(2, j);
      gik += stress(3) * (shapes.shapeDeriv(0, i) * shapes.shapeDeriv(1, j) +
                          shapes.shapeDeriv(1, i) * shapes.shapeDeriv(0, j));
      gik += stress(4) * (shapes.shapeDeriv(0, i) * shapes.shapeDeriv(2, j) +
                          shapes.shapeDeriv(2, i) * shapes.shapeDeriv(0, j));
      gik += stress(5) * (shapes.shapeDeriv(1, i) * shapes.shapeDeriv(2, j) +
                          shapes.shapeDeriv(2, i) * shapes.shapeDeriv(1, j));

      G.block(i * 3, j * 3, 3, 3) = gik * Types::Matrix33<prec>::Identity();
    }
  }

  return G;
}

auto EL303_ThermoMechanikSolid3D::getBTCBLinear(Geometry::H1Shapes &shapes,Types::MatrixXX<prec> &C, indexType nodeI, indexType nodeJ) -> Types::Matrix33<prec>
{

  Types::Matrix33<prec> res;

  prec N1x = shapes.shapeDeriv(0, nodeI);
  prec N1y = shapes.shapeDeriv(1, nodeI);
  prec N1z = shapes.shapeDeriv(2, nodeI);

  prec N2x = shapes.shapeDeriv(0, nodeJ);
  prec N2y = shapes.shapeDeriv(1, nodeJ);
  prec N2z = shapes.shapeDeriv(2, nodeJ);

  res(0,0) = N2x*(C(0,0)*N1x + C(0,3)*N1y + C(0,4)*N1z) + N2y*(C(0,3)*N1x + C(3,3)*N1y + C(3,4)*N1z) + N2z*(C(0,4)*N1x + C(3,4)*N1y + C(4,4)*N1z);
  res(0,1) = N2x*(C(0,3)*N1x + C(3,3)*N1y + C(3,4)*N1z) + N2y*(C(0,1)*N1x + C(1,3)*N1y + C(1,4)*N1z) + N2z*(C(0,5)*N1x + C(3,5)*N1y + C(4,5)*N1z);
  res(0,2) = N2x*(C(0,5)*N1x + C(3,5)*N1y + C(4,5)*N1z) + N2y*(C(0,4)*N1x + C(3,4)*N1y + C(4,4)*N1z) + N2z*(C(0,2)*N1x + C(2,3)*N1y + C(2,4)*N1z);
  res(1,0) = N2x*(C(0,1)*N1y + C(0,3)*N1x + C(0,5)*N1z) + N2y*(C(1,3)*N1y + C(3,3)*N1x + C(3,5)*N1z) + N2z*(C(1,4)*N1y + C(3,4)*N1x + C(4,5)*N1z);
  res(1,1) = N2x*(C(1,3)*N1y + C(3,3)*N1x + C(3,5)*N1z) + N2y*(C(1,1)*N1y + C(1,3)*N1x + C(1,5)*N1z) + N2z*(C(1,5)*N1y + C(3,5)*N1x + C(5,5)*N1z);
  res(1,2) = N2x*(C(1,5)*N1y + C(3,5)*N1x + C(5,5)*N1z) + N2y*(C(1,4)*N1y + C(3,4)*N1x + C(4,5)*N1z) + N2z*(C(1,2)*N1y + C(2,3)*N1x + C(2,5)*N1z);
  res(2,0) = N2x*(C(0,2)*N1z + C(0,4)*N1y + C(0,5)*N1x) + N2y*(C(2,3)*N1z + C(3,4)*N1y + C(3,5)*N1x) + N2z*(C(2,4)*N1z + C(4,4)*N1y + C(4,5)*N1x);
  res(2,1) = N2x*(C(2,3)*N1z + C(3,4)*N1y + C(3,5)*N1x) + N2y*(C(1,2)*N1z + C(1,4)*N1y + C(1,5)*N1x) + N2z*(C(2,5)*N1z + C(4,5)*N1y + C(5,5)*N1x);
  res(2,2) = N2x*(C(2,5)*N1z + C(4,5)*N1y + C(5,5)*N1x) + N2y*(C(2,4)*N1z + C(4,4)*N1y + C(4,5)*N1x) + N2z*(C(2,2)*N1z + C(2,4)*N1y + C(2,5)*N1x);

  return res;
}

void EL303_ThermoMechanikSolid3D::toParaviewAdaper(
    PointerCollection &pointers, FiniteElement::Volume &elem,
    vtkPlotInterface &paraviewAdapter, ParaviewSwitch control) {
  int matNum = static_cast<int>(elem.getMaterial()->getNumber());
  switch (control) {
  case ParaviewSwitch::Mesh: {
    elem.geometryToParaview(pointers, paraviewAdapter, 0, matNum);

  } break;
  case ParaviewSwitch::Solution: {
    elem.H1SolutionToParaview(pointers, paraviewAdapter, 0, matNum,
                               this->meshIdDisp, this->intOrderDisp,
                               paraviewNames::DisplacementName());
    elem.H1SolutionToParaview(pointers, paraviewAdapter, 0, matNum,
                               this->meshidTemp, this->intOrderDisp,
                               paraviewNames::TemperatureName());
  } break;
  case ParaviewSwitch::Weights: {
    elem.computeWeightsParaview(pointers, paraviewAdapter, 0, matNum);
  } break;
  case ParaviewSwitch::ProjectedValues: {
    auto GP = elem.getIntegrationPoints(pointers);
    GP.setOrder(this->intOrderDisp * 2);

    auto Hist = elem.getHistoryDataIterator(pointers);

    Materials::MaterialTransferData materialData;
    materialData.historyData = &Hist;
    materialData.strains.resize(6);
    std::vector<DegreeOfFreedom *> Dofs, TDofs;

    elem.getH1Dofs(pointers, Dofs, this->meshIdDisp, this->intOrderDisp);
    elem.getH1Dofs(pointers, TDofs, this->meshidTemp, this->intOrderDisp);
    auto solution = elem.getSolution(pointers, Dofs);
    auto solutionTemp = elem.getSolution(pointers, TDofs);
    indexType numDispDofs = Dofs.size();
    Types::VectorX<prec>  shapeMat;
    Types::Matrix3X<prec> tempGradient;
    shapeMat.resize(numDispDofs);
    shapeMat.setZero();
    tempGradient.resize(3, numDispDofs);
    tempGradient.setZero();

    for (auto i : GP) {
      auto jaco = elem.getJacobian(pointers, i);
      auto shapes = elem.getH1Shapes(pointers, this->intOrderDisp, jaco, i);
      for (auto j = 0; j < numDispDofs / 3; ++j) {
        shapeMat(3 * j) = shapes.shapes(j);
        tempGradient(0, 3 * j) = shapes.shapeDeriv(0, j);
        tempGradient(1, 3 * j) = shapes.shapeDeriv(1, j);
        tempGradient(2, 3 * j) = shapes.shapeDeriv(2, j);
      }

      if (this->mode == 1) {
        auto BMat = this->getLinearBMatrix(shapes);
        materialData.strains = BMat * solution;
      } else {
        auto F = this->getDeformationGradient(shapes, solution);
        auto E =
            (F.transpose() * F - Types::Matrix33<prec>::Identity()) * prec(0.5);
        materialData.strains(0) = E(0, 0);
        materialData.strains(1) = E(1, 1);
        materialData.strains(2) = E(2, 2);
        materialData.strains(3) = E(0, 1);
        materialData.strains(4) = E(0, 2);
        materialData.strains(5) = E(1, 2);
      }

      Types::VectorX<prec> TGrad = (tempGradient * solutionTemp).eval();
      Types::VectorX<prec> Flux = this->kappa * TGrad;

      elem.getMaterialFormulation(pointers)->getMaterialData(pointers, materialData,i);
      elem.projectDataToParaviewVertices(
          pointers, paraviewAdapter, 0, matNum, this->intOrderDisp, i,
          materialData.strains, 6, paraviewNames::strainName());
      elem.projectDataToParaviewVertices(
          pointers, paraviewAdapter, 0, matNum, this->intOrderDisp, i,
          materialData.stresses, 6, paraviewNames::stressName());

      
      elem.projectDataToParaviewVertices(
          pointers, paraviewAdapter, 0, matNum, this->intOrderDisp, i, TGrad, 3,
          paraviewNames::TemperatureGradientName());
      elem.projectDataToParaviewVertices(
          pointers, paraviewAdapter, 0, matNum, this->intOrderDisp, i, Flux, 3,
          paraviewNames::HeatFluxName());
      Hist.next();
    }
  }break;
  default:
    break;
  }
}

auto EL303_ThermoMechanikSolid3D::getHistoryDataStructure()
    -> const HistoryDataStructure & {
  return m_HistoryDataStructure;
}

auto EL303_ThermoMechanikSolid3D::getNumberOfIntergrationPoints(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> indexType {
  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(this->intOrderDisp * 2);
  return GP.getTotalGP();
}

const HistoryDataStructure EL303_ThermoMechanikSolid3D::m_HistoryDataStructure({},{});

} // namespace Elementformulations
} // namespace HierAMuS
