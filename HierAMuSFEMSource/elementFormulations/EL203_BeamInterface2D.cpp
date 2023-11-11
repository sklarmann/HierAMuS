// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include "geometry/Edges/EdgesData.h"
#include "geometry/GeometryTypes.h"
#include "geometry/VertexData.h"
#include "plot/vtkplotClassBase.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <iomanip>
#include <iostream>
#include <limits>
#include <ostream>


#include <elementFormulations/EL203_BeamInterface2D.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include "control/ParameterList.h"
#include <pointercollection/pointercollection.h>


#include <finiteElements/GenericFiniteElement.h>
#include <finiteElements/beamInterfaceElement2D.h>

#include <materials/GenericMaterialFormulation.h>
#include <materials/MaterialformulationList.h>

#include <geometry/GeometryBaseData.h>
#include "geometry/GeometryData.h"

#include <math/MatrixOperations.h>

#include <elementFormulations/GenericElementFormulation.h>
#include <solver/GenericSolutionState.h>

#include <types/MatrixTypes.h>

#include <Eigen/Eigenvalues>
#include <sstream>
#include <stdexcept>

#include <shapefunctions/LegendreShapes.h>
#include <vector>

#include <Timer.h>

#include <Eigen/Sparse>

#ifdef USE_VTK
#include <vtkCellType.h>
#endif

namespace HierAMuS::Elementformulations {

EL203_BeamInterface2D::EL203_BeamInterface2D(PointerCollection *ptrCol)
    : GenericElementFormulationInterface(ptrCol) {}

EL203_BeamInterface2D::~EL203_BeamInterface2D() {

  for (auto &i : this->Materials) {
    delete i;
  }
}

void EL203_BeamInterface2D::readData(PointerCollection &pointers,
                                     ParameterList &list) {
  this->meshIdDisp = list.getIndexVal("meshiddisp");
  this->meshIdWarp = list.getIndexVal("meshidwarp");
  this->meshIdRot = list.getIndexVal("meshidrot");
  this->mode = list.getIndexVal("mode");
  this->intOrder = list.getIndexVal("intorder");

  Types::MatrixXX<prec> tempMat, tempMat2;
  tempMat = list.getPrecMatrix("edgelist");
  tempMat2 = list.getPrecMatrix("matlist");

  if (tempMat.cols() == tempMat2.cols() && tempMat.rows() == tempMat2.rows()) {
    for (auto i = 0; i < tempMat.cols(); ++i) {
      for (auto j = 0; j < tempMat.rows(); ++j) {
        auto ednum = static_cast<indexType>(tempMat(j, i));
        auto manum = static_cast<indexType>(tempMat2(j, i));
        this->edgeMaterialMap[ednum] = manum;
      }
    }
  }

  this->beamVertexNumber = list.getIndexVal("vertex");

  if (this->materialList.size() != this->edgeList.size()) {
    throw std::runtime_error("Error in EL203, amount of material numbers does "
                             "not match amount of elements!");
  }
  
  auto &Log = pointers.getSPDLogger();

  Log.info("\n{:-<100}\n"
                "*   Element 201, specified Options\n"
                "    Mesh id for displacement nodes:         {:>12}\n"
                "    Shape function order for displacements: {:>12}\n"
                "{:-<100}\n",
                "",
                this->meshIdDisp,
                this->intOrder,
                "");


  this->messageUnprocessed(pointers, list, "EL203_BeamInterface2D");
}

void EL203_BeamInterface2D::setDegreesOfFreedom(
    PointerCollection &pointers, FiniteElement::beamInterfaceElement2D &elem) {

  elem.setMeshIdDisp(this->meshIdDisp);
  elem.setMeshIdRot(this->meshIdRot);
  elem.setMeshIdWarp(this->meshIdWarp);

  switch (this->mode) {
  case 1:
    this->setDegreesOfFreedomV1(pointers, &elem);
    break;
  case 2:
    this->setDegreesOfFreedomV2(pointers, &elem);
    break;
  case 3:
    this->setDegreesOfFreedomV3(pointers, &elem);
    break;
  case 4:
    this->setDegreesOfFreedomV4(pointers, &elem);
    break;
  case 5:
    this->setDegreesOfFreedomV5(pointers, elem);
    break;
  case 6:
    this->setDegreesOfFreedomV6(pointers, elem);
    break;
  case 7:
    this->setDegreesOfFreedomV6(pointers, elem);
    break;
  default:
    break;
  }

  // ee->setShapes(pointers, this->meshIdDisp, this->meshIdWarp,
  // this->meshIdRot, this->intOrder);
}

void EL203_BeamInterface2D::AdditionalOperations(
  PointerCollection &pointers, FiniteElement::beamInterfaceElement2D &elem) {
  
  switch (this->mode) {
  case 1:
    this->AdditionalOperationsV1(pointers, &elem);
    break;
  case 2:
    this->AdditionalOperationsV2(pointers, &elem);
    break;
  case 3:
    this->AdditionalOperationsV3(pointers, &elem);
    break;
  case 4:
    this->AdditionalOperationsV4(pointers, &elem);
    break;
  case 5:
    this->AdditionalOperationsV5(pointers, elem);
    break;
  case 6:
    this->AdditionalOperationsV6(pointers, elem);
    break;
  case 7:
    this->AdditionalOperationsV6(pointers, elem);
    break;
  default:
    break;
  }
}

auto EL203_BeamInterface2D::getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
-> std::vector<DegreeOfFreedom *> {
  std::vector<DegreeOfFreedom *> Dofs;

  auto *ee = dynamic_cast<FiniteElement::beamInterfaceElement2D *>(elem);

  std::vector<DegreeOfFreedom *> surfaceDispDofs;
  ee->getDofsOnSolid(pointers, surfaceDispDofs, this->meshIdDisp,
                       this->intOrder);
  std::vector<DegreeOfFreedom *> vertexDispDofs;
  ee->getDofsOnVert(pointers, vertexDispDofs, this->meshIdDisp);
  std::vector<DegreeOfFreedom *> vertexRotDofs;
  ee->getDofsOnVert(pointers, vertexRotDofs, this->meshIdRot);

  Dofs.clear();
  Dofs.insert(Dofs.end(), surfaceDispDofs.begin(), surfaceDispDofs.end());
  Dofs.insert(Dofs.end(), vertexDispDofs.begin(), vertexDispDofs.end());
  Dofs.insert(Dofs.end(), vertexRotDofs.begin(), vertexRotDofs.end());

  return Dofs;
}

void EL203_BeamInterface2D::setTangentResidual(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D &elem,
  Types::MatrixXX<prec> &stiffness, Types::VectorX<prec> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  auto *ee = dynamic_cast<FiniteElement::beamInterfaceElement2D *>(&elem);
  switch (this->mode) {
  case 1:
    this->setTangentResidualV1(pointers, ee, stiffness, residual, Dofs);
    break;
  case 2:
    this->setTangentResidualV2(pointers, ee, stiffness, residual, Dofs);
    break;
  case 3:
    this->setTangentResidualV3(pointers, ee, stiffness, residual, Dofs);
    break;
  case 4:
    this->setTangentResidualV4(pointers, ee, stiffness, residual, Dofs);
    break;
  case 5:
    this->setTangentResidualV5(pointers, ee, stiffness, residual, Dofs);
    break;
  case 6:
    this->setTangentResidualV6(pointers, ee, stiffness, residual, Dofs);
    break;
  case 7:
    this->setTangentResidualV7(pointers, ee, stiffness, residual, Dofs);
    break;
  default:
    break;
  }
}

void EL203_BeamInterface2D::setDegreesOfFreedomV1(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem) {
  elem->setDofsOnSolid(pointers, this->meshIdDisp, this->intOrder);
  elem->setDofsOnSolid(pointers, this->meshIdWarp, this->intOrder);
  elem->setDofsOnVert(pointers, this->meshIdDisp);
  elem->setDofsOnVert(pointers, this->meshIdRot);
}

void EL203_BeamInterface2D::AdditionalOperationsV1(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem) {
  // Geometry Set up
  elem->computeGeometry(pointers);

  // Degree of Freedom set up
  elem->setNodeMapDisp(pointers, this->meshIdDisp);
  elem->setNodeMapWarp(pointers, this->meshIdWarp);

  // Restrict unused nodes on Solid
  elem->setBCOnAllNodesSolid(pointers, this->meshIdDisp, 2);
  elem->setBCOnAllNodesSolid(pointers, this->meshIdWarp, 2);

  indexType numSolidNodes =
      elem->getNumberOfSolidNodes(pointers, this->meshIdDisp);
  elem->setBCOnNodeSolid(pointers, this->meshIdWarp, 0, 0);
  elem->setBCOnNodeSolid(pointers, this->meshIdWarp, 0, 1);
  elem->setBCOnNodeSolid(pointers, this->meshIdWarp, numSolidNodes - 1, 0);

  // Restrict degrees of freedom on vertex
  elem->setBCOnVert(pointers, this->meshIdDisp, 2);
  elem->setBCOnVert(pointers, this->meshIdRot, 1);
  elem->setBCOnVert(pointers, this->meshIdRot, 2);

  elem->computeShapesLocalWarping(pointers, this->meshIdWarp,
                                  this->intOrder);
  elem->computeSurfaceDispShapes(pointers, this->meshIdDisp,
                                 this->intOrder);
}

void EL203_BeamInterface2D::setTangentResidualV1(
  PointerCollection& pointers,
  FiniteElement::beamInterfaceElement2D *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {
  indexType numEdges = elem->getNumberOfEdges(pointers);

  Materials::MaterialTransferData materialData;

  IntegrationPoints GP = this->getIntegrationPoints(pointers, elem);

  Types::VectorX<prec> WarpingShapeX, WarpingShapeXY, WarpingShapeY,
      WarpingShapeYY, Nu, Nbeta, xiShape, dXiShape, surfShape, dSurfShape;

  Types::VectorX<prec> solution;

  elem->getSurfaceDispShapes(Nu, Nbeta);

  std::vector<DegreeOfFreedom *> surfaceDispDofs;
  elem->getDofsOnSolid(pointers, surfaceDispDofs, this->meshIdDisp,
                       this->intOrder);
  std::vector<DegreeOfFreedom *> surfaceWarpDofs;
  elem->getDofsOnSolid(pointers, surfaceWarpDofs, this->meshIdWarp,
                       this->intOrder);
  std::vector<DegreeOfFreedom *> vertexDispDofs;
  elem->getDofsOnVert(pointers, vertexDispDofs, this->meshIdDisp);
  std::vector<DegreeOfFreedom *> vertexRotDofs;
  elem->getDofsOnVert(pointers, vertexRotDofs, this->meshIdRot);

  Dofs.clear();
  Dofs.insert(Dofs.end(), surfaceDispDofs.begin(), surfaceDispDofs.end());
  Dofs.insert(Dofs.end(), surfaceWarpDofs.begin(), surfaceWarpDofs.end());
  Dofs.insert(Dofs.end(), vertexDispDofs.begin(), vertexDispDofs.end());
  Dofs.insert(Dofs.end(), vertexRotDofs.begin(), vertexRotDofs.end());

  elem->getSolution(pointers, Dofs, solution);

  stiffness.resize(Dofs.size(), Dofs.size());
  stiffness.setZero();
  residual.resize(Dofs.size());
  residual.setZero();

  Types::Matrix3X<prec> Bmat;
  Bmat.resize(3, Dofs.size());
  Bmat.setZero();

  for (auto i = 0; i < numEdges; ++i) {
    prec dA = elem->getJacXi() * elem->getJacEta(pointers, i);

    for (auto ngp = 0; ngp < GP.getTotalGP(); ++ngp) {
      GP.setCurrNumber(ngp);
      elem->getWarpingShapes(pointers, WarpingShapeX, WarpingShapeY,
                             WarpingShapeXY, WarpingShapeYY, i, GP.getEta(ngp),
                             this->intOrder, this->meshIdWarp);

      elem->getXiShape(pointers, this->intOrder, GP.getXi(ngp), xiShape,
                       dXiShape);
      elem->getLocalSurfaceDispShapesSorted(pointers, surfShape,
                                            dSurfShape, i, GP.getEta(ngp),
                                            this->intOrder, this->meshIdDisp);
      auto z = elem->getZCoordinate(pointers, i, GP.getEta(ngp));
      indexType pos = 0;
      for (auto j = 0; j < surfaceDispDofs.size() / 3; ++j) {
        Bmat(0, pos) = dXiShape(0) * (Nu(j) - z * Nbeta(j));
        Bmat(2, pos + 1) = -xiShape(0) * Nbeta(j) + dXiShape(0) * Nu(j);
        // Bmat(0,pos) = dXiShape(0)*surfShape(j);
        // Bmat(1,pos+1) = xiShape(0)*dSurfShape(j);
        // Bmat(2,pos) = xiShape(0)*dSurfShape(j);
        // Bmat(2,pos+1) = dXiShape(0)*surfShape(j);

        pos += 3;
      }

      for (auto j = 0; j < surfaceWarpDofs.size() / 3; ++j) {
        Bmat(1, pos + 1) = WarpingShapeYY(j);
        Bmat(2, pos) = WarpingShapeXY(j);
        pos += 3;
      }
      // std::cout << WarpingShapeXY << std::endl;
      // std::cout << WarpingShapeYY << std::endl;

      Bmat(0, pos) = dXiShape(1);
      Bmat(2, pos + 1) = dXiShape(1);
      pos += 3;
      Bmat(0, pos) = -z * dXiShape(1);
      Bmat(2, pos) = -xiShape(1);

      materialData.strains = Bmat * solution;
      elem->getMaterialFormulation(pointers)->getMaterialData(pointers, materialData, GP.getIntegrationPoint());

      // std::cout << GP.getXi(ngp) << " " << GP.getEta(ngp) <<
      // GP.getWeight(ngp) << std::endl; std::cout << Bmat << std::endl;

      stiffness += Bmat.transpose() * materialData.materialTangent * Bmat * dA *
                   GP.getWeight(ngp);
    }
  }
}

void EL203_BeamInterface2D::setDegreesOfFreedomV2(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem) {
  elem->setDofsOnSolid(pointers, this->meshIdDisp, this->intOrder);
  elem->setDofsOnSolid(pointers, this->meshIdWarp, this->intOrder);
  elem->setDofsOnVert(pointers, this->meshIdDisp);
  elem->setDofsOnVert(pointers, this->meshIdRot);

  elem->setMeshIdDisp(this->meshIdDisp);
  elem->setMeshIdRot(this->meshIdRot);
  elem->setMeshIdWarp(this->meshIdWarp);
}

void EL203_BeamInterface2D::AdditionalOperationsV2(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem) {
  // Geometry Set up
  elem->computeGeometry(pointers);

  // Degree of Freedom set up
  elem->setNodeMapDisp(pointers, this->meshIdDisp);
  elem->setNodeMapWarp(pointers, this->meshIdWarp);

  // Restrict unused nodes on Solid
  elem->setBCOnAllNodesSolid(pointers, this->meshIdDisp, 2);
  elem->setBCOnAllNodesSolid(pointers, this->meshIdWarp, 2);
  elem->setBCOnAllNodesSolid(pointers, this->meshIdWarp, 1);
  // ee->setBCOnAllNodesSolid(pointers, this->meshIdWarp, 0);
  indexType numSolidNodes =
      elem->getNumberOfSolidNodes(pointers, this->meshIdDisp);
  elem->setBCOnNodeSolid(pointers, this->meshIdWarp, 0, 0);
  elem->setBCOnNodeSolid(pointers, this->meshIdWarp, 0, 1);
  elem->setBCOnNodeSolid(pointers, this->meshIdWarp, 0, 2);
  elem->setBCOnNodeSolid(pointers, this->meshIdWarp, numSolidNodes - 1, 0);

  // Restrict degrees of freedom on vertex
  elem->setBCOnVert(pointers, this->meshIdDisp, 2);
  elem->setBCOnVert(pointers, this->meshIdRot, 1);
  elem->setBCOnVert(pointers, this->meshIdRot, 2);

  elem->computeShapesLocalWarping(pointers, this->meshIdWarp,
                                  this->intOrder);
  elem->computeSurfaceDispShapes(pointers, this->meshIdDisp,
                                 this->intOrder);
}

void EL203_BeamInterface2D::setTangentResidualV2(
  PointerCollection& pointers,
  FiniteElement::beamInterfaceElement2D *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {
  indexType numEdges = elem->getNumberOfEdges(pointers);

  Materials::MaterialTransferData materialData;

  IntegrationPoints GP = this->getIntegrationPoints(pointers,elem);

  Types::VectorX<prec> WarpingShapeX, WarpingShapeXY, WarpingShapeY,
      WarpingShapeYY, Nu, Nbeta, xiShape, dXiShape, surfShape, dSurfShape;
  Types::VectorX<prec> tempConst, tempLinear, solution;

  elem->getSurfaceDispShapes(Nu, Nbeta);

  // std::cout << Nu << std::endl;
  // std::cout << Nbeta << std::endl;

  std::vector<DegreeOfFreedom *> surfaceDispDofs;
  elem->getDofsOnSolid(pointers, surfaceDispDofs, this->meshIdDisp,
                       this->intOrder);
  std::vector<DegreeOfFreedom *> surfaceWarpDofs;
  elem->getDofsOnSolid(pointers, surfaceWarpDofs, this->meshIdWarp,
                       this->intOrder);
  std::vector<DegreeOfFreedom *> vertexDispDofs;
  elem->getDofsOnVert(pointers, vertexDispDofs, this->meshIdDisp);
  std::vector<DegreeOfFreedom *> vertexRotDofs;
  elem->getDofsOnVert(pointers, vertexRotDofs, this->meshIdRot);

  elem->getSolution(pointers, Dofs, solution);

  tempConst.resize(surfaceWarpDofs.size() / 3);
  tempLinear.resize(surfaceWarpDofs.size() / 3);
  tempConst.setZero();
  tempLinear.setZero();

  Dofs.clear();
  Dofs.insert(Dofs.end(), surfaceDispDofs.begin(), surfaceDispDofs.end());
  Dofs.insert(Dofs.end(), surfaceWarpDofs.begin(), surfaceWarpDofs.end());
  Dofs.insert(Dofs.end(), vertexDispDofs.begin(), vertexDispDofs.end());
  Dofs.insert(Dofs.end(), vertexRotDofs.begin(), vertexRotDofs.end());

  stiffness.resize(Dofs.size(), Dofs.size());
  stiffness.setZero();
  residual.resize(Dofs.size());
  residual.setZero();

  Types::Matrix3X<prec> Bmat;
  Bmat.resize(3, Dofs.size());
  Bmat.setZero();

  for (auto i = 0; i < numEdges; ++i) {
    prec dA = elem->getJacXi() * elem->getJacEta(pointers, i);
    for (auto ngp = 0; ngp < GP.getTotalGP(); ++ngp) {
      GP.setCurrNumber(ngp);
      Bmat.setZero();
      elem->getWarpingShapes(pointers, WarpingShapeX, WarpingShapeY,
                             WarpingShapeXY, WarpingShapeYY, i, GP.getEta(ngp),
                             this->intOrder, this->meshIdWarp);

      elem->getXiShape(pointers, this->intOrder, GP.getXi(ngp), xiShape,
                       dXiShape);
      elem->getLocalSurfaceDispShapesSorted(pointers, surfShape,
                                            dSurfShape, i, GP.getEta(ngp),
                                            this->intOrder, this->meshIdDisp);
      auto z = elem->getZCoordinate(pointers, i, GP.getEta(ngp));

      indexType pos = 0;
      for (auto j = 0; j < surfaceDispDofs.size() / 3; ++j) {

        // Shell kinematics
        Bmat(0, pos) = dXiShape(0) * surfShape(j);
        Bmat(1, pos + 1) = xiShape(0) * dSurfShape(j);
        Bmat(2, pos) = xiShape(0) * dSurfShape(j);
        Bmat(2, pos + 1) = dXiShape(0) * surfShape(j);

        // Beam kinematics
        // Bmat(0, pos) = dXiShape(0) * (Nu(j) - z * Nbeta(j)); // uxx
        // Bmat(2, pos) = -xiShape(0) * Nbeta(j); // uxy
        // Bmat(2, pos + 1) = dXiShape(0) * Nu(j); // uyx
        // Bmat(2,pos+1) = dXiShape(0)*Nu(j);

        pos += 3;
      }

      for (auto j = 0; j < surfaceWarpDofs.size() / 3; ++j) {

        // Shell kinematics
        Bmat(0, pos) = dXiShape(1) * WarpingShapeX(j);
        Bmat(1, pos + 1) = xiShape(1) * WarpingShapeYY(j);
        Bmat(2, pos) = xiShape(1) * WarpingShapeXY(j);
        Bmat(2, pos + 1) = dXiShape(1) * WarpingShapeY(j);

        // Constant Warping Function
        // Bmat(1, pos + 1) = WarpingShapeYY(j);
        // Bmat(2, pos) = WarpingShapeXY(j);
        pos += 3;
      }
      // std::cout << WarpingShapeXY << std::endl;
      // std::cout << WarpingShapeYY << std::endl;

      Bmat(0, pos) = dXiShape(1);
      Bmat(2, pos + 1) = dXiShape(1);
      pos += 3;
      Bmat(0, pos) = -z * dXiShape(1);
      Bmat(2, pos) = -xiShape(1);

      materialData.strains = Bmat * solution;

      elem->getMaterialFormulation(pointers)->getMaterialData(pointers, materialData,GP.getIntegrationPoint());

      // std::cout << GP.getXi(ngp) << " " << GP.getEta(ngp) <<
      // GP.getWeight(ngp) << std::endl; std::cout << Bmat << std::endl;
      stiffness += Bmat.transpose() * materialData.materialTangent * Bmat * dA *
                   GP.getWeight(ngp);
    }
  }
  // std::cout << tempConst << std::endl << std::endl;
  // std::cout << tempLinear << std::endl;
  // std::cout << stiffness.eigenvalues() << std::endl;
  //   elem->getWarpingShapes(pointers, WarpingShapeX, WarpingShapeY,
  //                        WarpingShapeXY, WarpingShapeYY, 1, 0,
  //                        this->meshIdWarp);

  // Types::VectorX<prec> shpx, shpy;
  // Types::Matrix2X<prec> dshpx, dshpy;
  // Types::MatrixXX<prec> cmat, cmatred;
  // Types::MatrixXX<prec> bmat, bmat2;

  // elem->getDofs(pointers, Dofs, this->meshIdDisp, this->meshIdWarp,
  // this->meshIdRot, this->intOrder); stiffness.resize(Dofs.size(),
  // Dofs.size()); residual.resize(Dofs.size()); stiffness.setZero();
  // residual.setZero();
  // bmat.resize(3, Dofs.size());
  // bmat.setZero();
  // bmat2.resize(3, Dofs.size());
  // bmat2.setZero();

  // cmat.resize(3,3);
  // cmatred.resize(3,3);
  // cmat.setZero();
  // cmatred.setZero();

  // if (this->planeStrain == 1) {
  //     prec fac =
  //         this->emod / ((prec)1 + this->nu) / ((prec)1 - (prec)2 * this->nu);
  //     cmat(0, 0) = ((prec)1 - this->nu);
  //     cmat(0, 1) = (this->nu);
  //     cmat(1, 1) = ((prec)1 - this->nu);
  //     cmat(1, 0) = (this->nu);
  //     cmat *= fac;
  //     cmat(2, 2) = this->emod / ((prec)1 + this->nu) / (prec)2;
  //     for(auto i=0;i<3;++i){
  //         cmatred(1,i) = cmat(1,i);
  //     }

  // }
  // else {
  //     prec fac = this->emod / ((prec)1 - this->nu * this->nu);
  //     cmat(0, 0) = ((prec)1);
  //     cmat(0, 1) = (this->nu);
  //     cmat(1, 1) = ((prec)1);
  //     cmat(1, 0) = (this->nu);
  //     cmat *= fac;
  //     cmat(2, 2) = this->emod / ((prec)1 + this->nu) / (prec)2;
  //     cmat *= this->thick;
  //     for(auto i=0;i<3;++i){
  //         cmatred(1,i) = cmat(1,i);
  //     }
  // }

  // std::vector<prec> xsi, eta, weight;
  // prec detj;
  // elem->getGaussPoints(3, weight, xsi, eta);
  // //linearGP(eta,weight,2);
  // //for (auto &i : weight)
  // //{
  // //    i *= prec(2);
  // //    xsi.push_back(prec(0));
  // //}
  // indexType nnodes = Dofs.size() / 3;
  // indexType normalnodes = nnodes-2;
  // normalnodes /= 2;
  // Types::VectorX<prec> etashp, detashp, xishp, dxishp;
  // for (auto edgeNum = 0; edgeNum < numEdges; ++edgeNum) {
  //     for (auto gp = 0; gp < xsi.size(); ++gp) {
  //         elem->getShapes(pointers, this->intOrder, this->meshIdDisp,
  //             this->meshIdWarp, this->meshIdRot,
  //             shpx, dshpx, shpy, dshpy, detj,
  //             xsi[gp], eta[gp], edgeNum);
  //         elem->getEtaShape(pointers, this->meshIdDisp, edgeNum,
  //         this->intOrder, eta[gp], etashp, detashp);
  //         elem->getXiShape(pointers, this->intOrder, xsi[gp], xishp,
  //         dxishp); prec y =
  //         elem->getZCoordinate(pointers,edgeNum,eta[gp]);
  //         //for (auto i = 0; i < normalnodes; ++i) {
  //         for (auto i = 0; i < nnodes; ++i) {
  //             bmat(0, 3 * i) = dshpx(0, i);
  //             bmat(1, 3 * i + 1) = dshpy(1, i);
  //             //bmat(1, 3 * i + 1) = detashp(i);

  //             bmat(2, 3 * i) = dshpx(1, i);
  //             bmat(2, 3 * i + 2) = dshpy(0, i);

  //             //uxy(3*i) = dshpx(1, i);
  //         }

  //         stiffness += bmat.transpose() * cmat * bmat * detj * weight[gp];
  //         //stiffness += uxy * lam.transpose() * detj * weight[gp];
  //         //stiffness += lam * uxy.transpose() * detj * weight[gp];
  //         //stiffness += prec(200000) * bmat2.transpose() * bmat2 * detj *
  //         weight[gp];
  //         //stiffness +=  prec(200) *
  //         bmat.block(1,0,1,Dofs.size()).transpose() * ((cmatred) *
  //         bmat).block(1,0,1,Dofs.size()) * detj * weight[gp];
  //         //stiffness += prec(10000000) * bmat.transpose() *
  //         cmatred.transpose() * cmatred * bmat * detj * weight[gp];
  //     }
  // }
  // std::cout << stiffness.eigenvalues() << std::endl;
}

void EL203_BeamInterface2D::setDegreesOfFreedomV3(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem) {
  elem->setDofsOnSolid(pointers, this->meshIdDisp, this->intOrder);
  elem->setDofsOnSolid(pointers, this->meshIdWarp, this->intOrder);
  elem->setDofsOnVert(pointers, this->meshIdDisp);
  elem->setDofsOnVert(pointers, this->meshIdRot);

  elem->setMeshIdDisp(this->meshIdDisp);
  elem->setMeshIdRot(this->meshIdRot);
  elem->setMeshIdWarp(this->meshIdWarp);
}

void EL203_BeamInterface2D::AdditionalOperationsV3(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem) {
  // Geometry Set up
  elem->computeGeometry(pointers);

  // Degree of Freedom set up
  elem->setNodeMapDisp(pointers, this->meshIdDisp);
  elem->setNodeMapWarp(pointers, this->meshIdWarp);

  // Restrict unused nodes on Solid
  elem->setBCOnAllNodesSolid(pointers, this->meshIdDisp, 2);
  elem->setBCOnAllNodesSolid(pointers, this->meshIdWarp, 2);
  elem->setBCOnAllNodesSolid(pointers, this->meshIdWarp, 1);
  elem->setBCOnAllNodesSolid(pointers, this->meshIdWarp, 0);

  // Restrict degrees of freedom on vertex
  elem->setBCOnVert(pointers, this->meshIdDisp, 2);
  elem->setBCOnVert(pointers, this->meshIdRot, 1);
  elem->setBCOnVert(pointers, this->meshIdRot, 2);

  elem->computeShapesLocalWarping(pointers, this->meshIdWarp,
                                  this->intOrder);
  elem->computeSurfaceDispShapes(pointers, this->meshIdDisp,
                                 this->intOrder);
}

void EL203_BeamInterface2D::setTangentResidualV3(
  PointerCollection& pointers,
  FiniteElement::beamInterfaceElement2D *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {
  indexType numEdges = elem->getNumberOfEdges(pointers);

  Materials::MaterialTransferData materialData;

  IntegrationPoints GP = this->getIntegrationPoints(pointers, elem);

  Types::VectorX<prec> WarpingShapeX, WarpingShapeXY, WarpingShapeY,
      WarpingShapeYY, Nu, Nbeta, xiShape, dXiShape, surfShape, dSurfShape;
  Types::VectorX<prec> tempConst, tempLinear, solution;

  elem->getSurfaceDispShapes(Nu, Nbeta);

  // std::cout << Nu << std::endl;
  // std::cout << Nbeta << std::endl;

  std::vector<DegreeOfFreedom *> surfaceDispDofs;
  elem->getDofsOnSolid(pointers, surfaceDispDofs, this->meshIdDisp,
                       this->intOrder);
  std::vector<DegreeOfFreedom *> vertexDispDofs;
  elem->getDofsOnVert(pointers, vertexDispDofs, this->meshIdDisp);
  std::vector<DegreeOfFreedom *> vertexRotDofs;
  elem->getDofsOnVert(pointers, vertexRotDofs, this->meshIdRot);

  tempConst.resize(surfaceDispDofs.size() / 3);
  tempLinear.resize(surfaceDispDofs.size() / 3);
  tempConst.setZero();
  tempLinear.setZero();

  Dofs.clear();
  Dofs.insert(Dofs.end(), surfaceDispDofs.begin(), surfaceDispDofs.end());
  Dofs.insert(Dofs.end(), vertexDispDofs.begin(), vertexDispDofs.end());
  Dofs.insert(Dofs.end(), vertexRotDofs.begin(), vertexRotDofs.end());

  elem->getSolution(pointers, Dofs, solution);

  stiffness.resize(Dofs.size(), Dofs.size());
  stiffness.setZero();
  residual.resize(Dofs.size());
  residual.setZero();

  indexType numCols, numRows;

  Types::Matrix3X<prec> Bmat, Amat, Amat2;
  numRows = 3;
  numCols = Dofs.size();
  Bmat.resize(numRows, numCols);
  Bmat.setZero();

  numCols = 3 + surfaceDispDofs.size() - 5;
  Amat.resize(numRows, numCols);
  Amat.setZero();
  numCols -= 3;
  Amat2.resize(numRows, numCols);
  Amat2.setZero();

  Types::MatrixXX<prec> Gmat, Hmat, Lmat, Jmat;
  Gmat.resize(Bmat.cols(), Amat.cols());
  Gmat.setZero();
  Hmat.resize(Amat.cols(), Amat.cols());
  Hmat.setZero();
  Lmat.resize(Amat.cols(), Amat2.cols());
  Lmat.setZero();
  Jmat.resize(Amat.cols(), Amat.cols());
  Jmat.setZero();

  for (auto i = 0; i < numEdges; ++i) {
    prec dA = elem->getJacXi() * elem->getJacEta(pointers, i);
    for (auto ngp = 0; ngp < GP.getTotalGP(); ++ngp) {
      GP.setCurrNumber(ngp);
      Bmat.setZero();
      elem->getWarpingShapes(pointers, WarpingShapeX, WarpingShapeY,
                             WarpingShapeXY, WarpingShapeYY, i, GP.getEta(ngp),
                             this->intOrder, this->meshIdWarp);

      elem->getXiShape(pointers, this->intOrder, GP.getXi(ngp), xiShape,
                       dXiShape);
      elem->getLocalSurfaceDispShapesSorted(pointers, surfShape,
                                            dSurfShape, i, GP.getEta(ngp),
                                            this->intOrder, this->meshIdDisp);
      auto z = elem->getZCoordinate(pointers, i, GP.getEta(ngp));

      indexType pos = 0;
      for (auto j = 0; j < surfaceDispDofs.size() / 3; ++j) {
        Bmat(0, pos) = dXiShape(0) * surfShape(j);
        Bmat(1, pos + 1) = xiShape(0) * dSurfShape(j);
        Bmat(2, pos) = xiShape(0) * dSurfShape(j);
        Bmat(2, pos + 1) = dXiShape(0) * surfShape(j);

        pos += 3;
      }
      Bmat(0, pos) = dXiShape(1);
      Bmat(2, pos + 1) = dXiShape(1);
      pos += 3;
      Bmat(0, pos) = -z * dXiShape(1);
      Bmat(2, pos) = -xiShape(1);

      pos = 0;
      for (auto j = 0; j < surfaceDispDofs.size() / 3 - 2; ++j) {
        Amat2(0, pos) = WarpingShapeX(j);
        ++pos;
      }
      for (auto j = 0; j < surfaceDispDofs.size() / 3 - 1; ++j) {
        Amat2(1, pos) = WarpingShapeYY(j);
        ++pos;
      }
      for (auto j = 0; j < surfaceDispDofs.size() / 3 - 2; ++j) {
        Amat2(2, pos) = WarpingShapeXY(j);
        ++pos;
      }

      materialData.strains = Bmat * solution;

      elem->getMaterialFormulation(pointers)->getMaterialData(pointers, materialData,GP.getIntegrationPoint());

      Amat.block(0, 3, 3, Amat2.cols()) = Amat2;
      Amat(0, 0) = 1;
      Amat(0, 1) = z;
      Amat(2, 2) = 1;
      Gmat += Bmat.transpose() * Amat * dA * GP.getWeight(ngp);
      Jmat += Amat.transpose() * materialData.materialTangent * Amat * dA *
              GP.getWeight(ngp);
      Hmat += Amat.transpose() * Amat * dA * GP.getWeight(ngp);

      Lmat += Amat.transpose() * Amat2 * dA * GP.getWeight(ngp);
    }
  }
  auto Mmat = Hmat * Jmat.inverse() * Hmat.transpose();
  auto Mmatinv = Mmat.inverse();
  stiffness = Gmat * Mmatinv * Gmat.transpose() -
              Gmat * Mmatinv * Lmat *
                  (Lmat.transpose() * Mmat.inverse() * Lmat).inverse() *
                  Lmat.transpose() * Mmatinv * Gmat.transpose();
  // Eigen::EigenSolver<Types::MatrixXX<prec>> es(stiffness);
  // std::cout << "Eigenvalues of J are: " << es.eigenvalues() << std::endl;
}

void EL203_BeamInterface2D::setDegreesOfFreedomV4(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem) {
  elem->setDofsOnSolid(pointers, this->meshIdDisp, this->intOrder);
  elem->setDofsOnSolid(pointers, this->meshIdWarp, this->intOrder);
  elem->setDofsOnVert(pointers, this->meshIdDisp);
  elem->setDofsOnVert(pointers, this->meshIdRot);

  elem->setMeshIdDisp(this->meshIdDisp);
  elem->setMeshIdRot(this->meshIdRot);
  elem->setMeshIdWarp(this->meshIdWarp);
}

void EL203_BeamInterface2D::AdditionalOperationsV4(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem) {
  // Geometry Set up
  elem->computeGeometry(pointers);

  // Degree of Freedom set up
  elem->setNodeMapDisp(pointers, this->meshIdDisp);
  elem->setNodeMapWarp(pointers, this->meshIdWarp);

  // Restrict unused nodes on Solid
  elem->setBCOnAllNodesSolid(pointers, this->meshIdDisp, 2);
  elem->setBCOnAllNodesSolid(pointers, this->meshIdWarp, 2);
  elem->setBCOnAllNodesSolid(pointers, this->meshIdWarp, 1);
  elem->setBCOnAllNodesSolid(pointers, this->meshIdWarp, 0);

  // Restrict degrees of freedom on vertex
  elem->setBCOnVert(pointers, this->meshIdDisp, 2);
  elem->setBCOnVert(pointers, this->meshIdRot, 1);
  elem->setBCOnVert(pointers, this->meshIdRot, 2);

  elem->computeShapesLocalWarping(pointers, this->meshIdWarp,
                                  this->intOrder);
  elem->computeSurfaceDispShapes(pointers, this->meshIdDisp,
                                 this->intOrder);
}

void EL203_BeamInterface2D::setTangentResidualV4(
  PointerCollection& pointers,
  FiniteElement::beamInterfaceElement2D *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {
  indexType numEdges = elem->getNumberOfEdges(pointers);

  Materials::MaterialTransferData materialData;

  IntegrationPoints GP = this->getIntegrationPoints(pointers, elem);

  Types::VectorX<prec> WarpingShapeX, WarpingShapeXY, WarpingShapeY,
      WarpingShapeYY, Nu, Nbeta, xiShape, dXiShape, surfShape, dSurfShape;
  Types::VectorX<prec> tempConst, tempLinear, solution;

  elem->getSurfaceDispShapes(Nu, Nbeta);

  // std::cout << Nu << std::endl;
  // std::cout << Nbeta << std::endl;

  std::vector<DegreeOfFreedom *> surfaceDispDofs;
  elem->getDofsOnSolid(pointers, surfaceDispDofs, this->meshIdDisp,
                       this->intOrder);
  std::vector<DegreeOfFreedom *> vertexDispDofs;
  elem->getDofsOnVert(pointers, vertexDispDofs, this->meshIdDisp);
  std::vector<DegreeOfFreedom *> vertexRotDofs;
  elem->getDofsOnVert(pointers, vertexRotDofs, this->meshIdRot);

  elem->getSolution(pointers, Dofs, solution);

  tempConst.resize(surfaceDispDofs.size() / 3);
  tempLinear.resize(surfaceDispDofs.size() / 3);
  tempConst.setZero();
  tempLinear.setZero();

  Dofs.clear();
  Dofs.insert(Dofs.end(), surfaceDispDofs.begin(), surfaceDispDofs.end());
  Dofs.insert(Dofs.end(), vertexDispDofs.begin(), vertexDispDofs.end());
  Dofs.insert(Dofs.end(), vertexRotDofs.begin(), vertexRotDofs.end());

  stiffness.resize(Dofs.size(), Dofs.size());
  stiffness.setZero();
  residual.resize(Dofs.size());
  residual.setZero();

  indexType numCols, numRows;

  Types::Matrix3X<prec> Bmat, Amat, Amat2, AmatLocal;
  numRows = 3;
  numCols = Dofs.size();
  Bmat.resize(numRows, numCols);
  Bmat.setZero();

  numCols = 3 + surfaceDispDofs.size() - 5;
  Amat.resize(numRows, numCols);
  Amat.setZero();
  numCols -= 3;
  Amat2.resize(numRows, numCols);
  Amat2.setZero();

  AmatLocal.resize(3, numEdges * 3);
  AmatLocal.setZero();

  Types::MatrixXX<prec> Gmat, Hmat, Lmat, Jmat;
  Gmat.resize(Bmat.cols(), AmatLocal.cols());
  Gmat.setZero();
  Hmat.resize(AmatLocal.cols(), AmatLocal.cols());
  Hmat.setZero();
  Lmat.resize(AmatLocal.cols(), Amat2.cols());
  Lmat.setZero();
  Jmat.resize(AmatLocal.cols(), AmatLocal.cols());
  Jmat.setZero();

  for (auto i = 0; i < numEdges; ++i) {
    prec dA = elem->getJacXi() * elem->getJacEta(pointers, i);
    for (auto ngp = 0; ngp < GP.getTotalGP(); ++ngp) {
      GP.setCurrNumber(ngp);
      Bmat.setZero();
      elem->getWarpingShapes(pointers, WarpingShapeX, WarpingShapeY,
                             WarpingShapeXY, WarpingShapeYY, i, GP.getEta(ngp),
                             this->intOrder, this->meshIdWarp);

      elem->getXiShape(pointers, this->intOrder, GP.getXi(ngp), xiShape,
                       dXiShape);
      elem->getLocalSurfaceDispShapesSorted(pointers, surfShape,
                                            dSurfShape, i, GP.getEta(ngp),
                                            this->intOrder, this->meshIdDisp);
      auto z = elem->getZCoordinate(pointers, i, GP.getEta(ngp));

      indexType pos = 0;
      for (auto j = 0; j < surfaceDispDofs.size() / 3; ++j) {
        Bmat(0, pos) = dXiShape(0) * surfShape(j);
        Bmat(1, pos + 1) = xiShape(0) * dSurfShape(j);
        Bmat(2, pos) = xiShape(0) * dSurfShape(j);
        Bmat(2, pos + 1) = dXiShape(0) * surfShape(j);

        pos += 3;
      }
      Bmat(0, pos) = dXiShape(1);
      Bmat(2, pos + 1) = dXiShape(1);
      pos += 3;
      Bmat(0, pos) = -z * dXiShape(1);
      Bmat(2, pos) = -xiShape(1);

      pos = 0;
      for (auto j = 0; j < surfaceDispDofs.size() / 3 - 2; ++j) {
        Amat2(0, pos) = WarpingShapeYY(j);
        ++pos;
      }
      for (auto j = 0; j < surfaceDispDofs.size() / 3 - 1; ++j) {
        Amat2(1, pos) = WarpingShapeYY(j);
        ++pos;
      }
      for (auto j = 0; j < surfaceDispDofs.size() / 3 - 2; ++j) {
        Amat2(2, pos) = WarpingShapeXY(j);
        ++pos;
      }

      AmatLocal.setZero();

      for (auto nn = 0; nn < 3; ++nn) {
        AmatLocal(nn, 3 * i + nn) = prec(1);
      }

      elem->getMaterialFormulation(pointers)->getMaterialData(pointers, materialData,GP.getIntegrationPoint());

      Amat.block(0, 3, 3, Amat2.cols()) = Amat2;
      Amat(0, 0) = prec(1);
      Amat(0, 1) = z;
      Amat(2, 2) = prec(1);
      Gmat += Bmat.transpose() * AmatLocal * dA * GP.getWeight(ngp);
      Jmat += AmatLocal.transpose() * materialData.materialTangent * AmatLocal *
              dA * GP.getWeight(ngp);
      Hmat += AmatLocal.transpose() * AmatLocal * dA * GP.getWeight(ngp);
      Lmat += AmatLocal.transpose() * Amat2 * dA * GP.getWeight(ngp);
    }
  }
  auto Mmat = Hmat * Jmat.inverse() * Hmat.transpose();
  auto Mmatinv = Mmat.inverse();
  stiffness = Gmat * Mmatinv * Gmat.transpose() -
              Gmat * Mmatinv * Lmat *
                  (Lmat.transpose() * Mmat.inverse() * Lmat).inverse() *
                  Lmat.transpose() * Mmatinv * Gmat.transpose();
  // stiffness = Gmat*Mmatinv*Gmat.transpose() - Gmat*Mmatinv*Gmat.transpose();
  // Eigen::EigenSolver<Types::MatrixXX<prec>> es(stiffness);
  // std::cout << "Eigenvalues of J are: " << es.eigenvalues() << std::endl;
}

void EL203_BeamInterface2D::setDegreesOfFreedomV5(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D &elem) {
  elem.setDofsOnSolid(pointers, this->meshIdDisp, this->intOrder);
  elem.setDofsOnVert(pointers, this->meshIdDisp);
  elem.setDofsOnVert(pointers, this->meshIdRot);
}

void EL203_BeamInterface2D::AdditionalOperationsV5(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D &elem) {
  elem.setBCOnAllNodesSolid(pointers, this->meshIdDisp, 2);

  // Restrict degrees of freedom on vertex
  elem.setBCOnVert(pointers, this->meshIdDisp, 2);
  elem.setBCOnVert(pointers, this->meshIdRot, 1);
  elem.setBCOnVert(pointers, this->meshIdRot, 2);

  elem.computeWarpingShapesNew(pointers, this->intOrder);
}

void EL203_BeamInterface2D::setTangentResidualV5(
  PointerCollection& pointers,
  FiniteElement::beamInterfaceElement2D *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {
  auto GP = this->getIntegrationPoints(pointers, elem);
  indexType numberOfEdges = this->edgeMaterialMap.size();

  Types::VectorX<prec> WarpingShapeX, WarpingShapeXY, WarpingShapeY;
  Types::VectorX<prec> solution, xiShape, dXiShape;
  Types::VectorX<prec> surfShape, dSurfShape;

  std::vector<DegreeOfFreedom *> surfaceDispDofs;
  elem->getDofsOnSolid(pointers, surfaceDispDofs, this->meshIdDisp,
                       this->intOrder);
  std::vector<DegreeOfFreedom *> vertexDispDofs;
  elem->getDofsOnVert(pointers, vertexDispDofs, this->meshIdDisp);
  std::vector<DegreeOfFreedom *> vertexRotDofs;
  elem->getDofsOnVert(pointers, vertexRotDofs, this->meshIdRot);

  Dofs.clear();
  Dofs.insert(Dofs.end(), surfaceDispDofs.begin(), surfaceDispDofs.end());
  Dofs.insert(Dofs.end(), vertexDispDofs.begin(), vertexDispDofs.end());
  Dofs.insert(Dofs.end(), vertexRotDofs.begin(), vertexRotDofs.end());

  elem->getSolution(pointers, Dofs, solution);

  indexType numDofs;
  numDofs = Dofs.size();
  stiffness.resize(numDofs, numDofs);
  residual.resize(numDofs);
  stiffness.setZero();
  residual.setZero();

  indexType numCols, numRows;

  Types::Matrix3X<prec> Bmat, Amat, Amat2, AmatStress;
  numRows = 3;
  numCols = Dofs.size();
  Bmat.resize(numRows, numCols);
  Bmat.setZero();

  numCols = 3 + surfaceDispDofs.size() - 5;
  Amat.resize(numRows, numCols);
  Amat.setZero();
  AmatStress.resize(numRows, numCols - 1);
  AmatStress.setZero();
  numCols -= 3;
  Amat2.resize(numRows, numCols);
  Amat2.setZero();

  Types::MatrixXX<prec> Gmat, Hmat, Lmat, Jmat;
  Gmat.resize(Bmat.cols(), AmatStress.cols());
  Gmat.setZero();
  Hmat.resize(AmatStress.cols(), Amat.cols());
  Hmat.setZero();
  Lmat.resize(AmatStress.cols(), Amat2.cols());
  Lmat.setZero();
  Jmat.resize(Amat.cols(), Amat.cols());
  Jmat.setZero();

  Types::Vector3<prec> A1;
  A1 = elem->getA1(pointers);
  Types::Vector3<prec> A2;
  A2 = elem->getA2(pointers);

  for (auto edgeNum = 0; edgeNum < numberOfEdges; ++edgeNum) {
    indexType gEdge = elem->getGlobalEdgeNumber(pointers, edgeNum);
    indexType mNum = this->edgeMaterialMap[gEdge];
    auto material = pointers.getMaterialFormulationList()->getMaterial(mNum);

    Materials::MaterialTransferData matData;

    for (auto ngp = 0; ngp < GP.getTotalGP(); ++ngp) {
      GP.setCurrNumber(ngp);
      prec dA =
          elem->getDA(pointers, edgeNum, GP.getXi(ngp), GP.getEta(ngp));
      Bmat.setZero();
      elem->getLocalWarpingShapesA1(pointers, WarpingShapeX,
                                    WarpingShapeXY, edgeNum, GP.getEta(ngp));
      elem->getXiShape(pointers, 1, GP.getXi(ngp), xiShape, dXiShape);
      elem->getLocalSurfaceDispShapesSorted(pointers, surfShape,
                                            dSurfShape, edgeNum, GP.getEta(ngp),
                                            this->intOrder, this->meshIdDisp);
      prec z = elem->getZCoordinate(pointers, edgeNum, GP.getEta(ngp));

      indexType pos = 0;
      for (auto j = 0; j < surfaceDispDofs.size() / 3; ++j) {
        Bmat.block(0, pos, 1, 3) = A1.transpose() * dXiShape(0) * surfShape(j);
        Bmat.block(1, pos, 1, 3) = A2.transpose() * xiShape(0) * dSurfShape(j);
        Bmat.block(2, pos, 1, 3) = A1.transpose() * xiShape(0) * dSurfShape(j);
        Bmat.block(2, pos, 1, 3) += A2.transpose() * dXiShape(0) * surfShape(j);
        // Bmat(0, pos) = dXiShape(0) * surfShape(j);
        // Bmat(1, pos + 1) = xiShape(0) * dSurfShape(j);
        // Bmat(2, pos) = xiShape(0) * dSurfShape(j);
        // Bmat(2, pos + 1) = dXiShape(0) * surfShape(j);

        pos += 3;
      }
      Bmat.block(0, pos, 1, 3) = A1.transpose() * dXiShape(1);
      Bmat.block(2, pos, 1, 3) = A2.transpose() * dXiShape(1);
      // Bmat(0, pos) = dXiShape(1);
      // Bmat(2, pos + 1) = dXiShape(1);
      pos += 3;
      Bmat(0, pos) = -z * dXiShape(1);
      Bmat(2, pos) = -xiShape(1);

      pos = 0;
      AmatStress.setZero();
      AmatStress(1, pos + edgeNum) = prec(1);
      for (auto j = 0; j < surfaceDispDofs.size() / 3 - 1; ++j) {
        // AmatStress(0, pos) = surfShape(j);
        ++pos;
      }

      AmatStress(1, pos + edgeNum) = prec(1);
      for (auto j = 0; j < surfaceDispDofs.size() / 3 - 1; ++j) {
        ++pos;
      }
      AmatStress(2, pos + edgeNum) = prec(1);
      for (auto j = 0; j < surfaceDispDofs.size() / 3 - 1; ++j) {
        ++pos;
      }

      pos = 0;
      for (auto j = 0; j < surfaceDispDofs.size() / 3 - 2; ++j) {
        Amat2(0, pos) = WarpingShapeX(j);
        ++pos;
      }
      for (auto j = 0; j < surfaceDispDofs.size() / 3 - 1; ++j) {
        Amat2(1, pos) = dSurfShape(j);
        ++pos;
      }
      for (auto j = 0; j < surfaceDispDofs.size() / 3 - 2; ++j) {
        Amat2(2, pos) = WarpingShapeXY(j);
        ++pos;
      }

      matData.strains = Bmat * solution;
      material->getMaterialData(pointers, matData,GP.getIntegrationPoint());

      dA *= prec(2);

      Amat.block(0, 3, 3, Amat2.cols()) = Amat2;
      Amat(0, 0) = 1;
      Amat(0, 1) = z;
      Amat(2, 2) = 1;
      Gmat += Bmat.transpose() * AmatStress * dA * GP.getWeight(ngp);
      Jmat += Amat.transpose() * matData.materialTangent * Amat * dA *
              GP.getWeight(ngp);
      Hmat += AmatStress.transpose() * Amat * dA * GP.getWeight(ngp);

      Lmat += AmatStress.transpose() * Amat2 * dA * GP.getWeight(ngp);
    }
  }
  auto Mmat = Hmat * Jmat.inverse() * Hmat.transpose();
  auto Mmatinv = Mmat.inverse();
  stiffness = Gmat * Mmatinv * Gmat.transpose() -
              Gmat * Mmatinv * Lmat *
                  (Lmat.transpose() * Mmatinv * Lmat).inverse() *
                  Lmat.transpose() * Mmatinv * Gmat.transpose();

  stiffness = (stiffness + stiffness.transpose()) / prec(2);
  // Eigen::EigenSolver<Types::MatrixXX<prec>> es(stiffness);
  // std::cout << "Eigenvalues of J are: " << es.eigenvalues() << std::endl;
  residual = stiffness * solution;
}

void EL203_BeamInterface2D::setDegreesOfFreedomV6(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D &elem) {
  elem.setDofsOnSolid(pointers, this->meshIdDisp, this->intOrder);
  elem.setDofsOnVert(pointers, this->meshIdDisp);
  elem.setDofsOnVert(pointers, this->meshIdRot);
}

void EL203_BeamInterface2D::AdditionalOperationsV6(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D &elem) {
  elem.setBCOnAllNodesSolid(pointers, this->meshIdDisp, 2);

  // Restrict degrees of freedom on vertex
  elem.setBCOnVert(pointers, this->meshIdDisp, 2);
  elem.setBCOnVert(pointers, this->meshIdRot, 1);
  elem.setBCOnVert(pointers, this->meshIdRot, 2);

  elem.computeWarpingShapesNew(pointers, this->intOrder);

  // Calculation data which need to be stored
  indexType numberOfEdges = this->edgeMaterialMap.size();
  indexType numCols = (this->intOrder + 1) + this->intOrder * 2;
  numCols *= numberOfEdges;
  elem.request_element_data_field(pointers, m_resEps_id, numCols, 1);
  elem.request_element_data_field(pointers, m_resSig_id, numCols, 1);
  elem.request_element_data_field(pointers, m_Eps_id, numCols, 1);
  elem.request_element_data_field(pointers, m_Sig_id, numCols, 1);

  indexType numCols2 = numberOfEdges * 3 * this->intOrder - 2;
  elem.request_element_data_field(pointers, m_resEpsC_id, numCols2, 1);
  elem.request_element_data_field(pointers, m_EpsC_id, numCols2, 1);


  std::vector<DegreeOfFreedom *> surfaceDispDofs;
  elem.getDofsOnSolid(pointers, surfaceDispDofs, this->meshIdDisp,
                       this->intOrder);
  indexType dofCount = surfaceDispDofs.size() + 6;
  elem.request_element_data_field(pointers,m_G_id,dofCount,numCols);
  elem.request_element_data_field(pointers, m_H_id, numCols, numCols);
  elem.request_element_data_field(pointers, m_J_id, numCols, numCols);
  elem.request_element_data_field(pointers, m_L_id, numCols, numCols2);
  elem.request_element_data_field(pointers, m_M_id, numCols, numCols);
  
  


  // Compute matrices
  Types::MatrixXX<prec> &Gmat = elem.get_element_data_field(pointers, m_G_id);
  Types::MatrixXX<prec> &Hmat = elem.get_element_data_field(pointers, m_H_id);
  Types::MatrixXX<prec> &Lmat = elem.get_element_data_field(pointers, m_L_id);
  Types::MatrixXX<prec> &Jmat = elem.get_element_data_field(pointers, m_J_id);

  Types::VectorX<prec> WarpingShapeX, WarpingShapeXY, WarpingShapeY;
  Types::VectorX<prec> xiShape, dXiShape;
  Types::VectorX<prec> surfShape, dSurfShape;

  
  Types::Matrix3X<prec> AmatStress =
      Types::Matrix3X<prec>::Zero(3, numCols);

  std::vector<Eigen::Triplet<prec, indexType>> triplets;
  triplets.reserve((this->intOrder) * 3);
  Eigen::SparseMatrix<prec> AmatSparseSig(AmatStress.rows(), AmatStress.cols());

  Types::Matrix3X<prec> AmatEpsC =
      Types::Matrix3X<prec>::Zero(3, numCols2);

  auto GP = this->getIntegrationPoints(pointers,&elem);
  auto History = elem.getHistoryDataIterator(pointers);
  Types::Matrix3X<prec> Bmat(3, dofCount);

  
  Types::Vector3<prec> A1 = elem.getA1(pointers);
  Types::Vector3<prec> A2 = elem.getA2(pointers);
  Materials::MaterialTransferData matData;
  matData.historyData = &History;
  matData.strains.resize(3);

  for (auto edgeNum = 0; edgeNum < numberOfEdges; ++edgeNum) {
    indexType gEdge = elem.getGlobalEdgeNumber(pointers, edgeNum);
    indexType mNum = this->edgeMaterialMap[gEdge];
    auto material = pointers.getMaterialFormulationList()->getMaterial(mNum);

    

    for (auto ngp = 0; ngp < GP.getTotalGP(); ++ngp) {
      GP.setCurrNumber(ngp);
      prec dA = elem.getDA(pointers, edgeNum, GP.getXi(ngp), GP.getEta(ngp));
      Bmat.setZero();
      elem.getLocalWarpingShapesA1(pointers, WarpingShapeX, WarpingShapeXY,
                                    edgeNum, GP.getEta(ngp));
      elem.getXiShape(pointers, this->intOrder, GP.getXi(ngp), xiShape,
                       dXiShape);
      elem.getLocalSurfaceDispShapesSorted(pointers, surfShape, dSurfShape,
                                            edgeNum, GP.getEta(ngp),
                                            this->intOrder, this->meshIdDisp);
      prec z = elem.getZCoordinate(pointers, edgeNum, GP.getEta(ngp));

      indexType pos = 0;
      for (auto j = 0; j < surfaceDispDofs.size() / 3; ++j) {
        Bmat.block(0, pos, 1, 3) = A1.transpose() * dXiShape(0) * surfShape(j);
        Bmat.block(1, pos, 1, 3) = A2.transpose() * xiShape(0) * dSurfShape(j);
        Bmat.block(2, pos, 1, 3) = A1.transpose() * xiShape(0) * dSurfShape(j);
        Bmat.block(2, pos, 1, 3) += A2.transpose() * dXiShape(0) * surfShape(j);

        pos += 3;
      }
      Bmat.block(0, pos, 1, 3) = A1.transpose() * dXiShape(1);
      Bmat.block(2, pos, 1, 3) = A2.transpose() * dXiShape(1);
      pos += 3;
      Bmat(0, pos) = -z * dXiShape(1);
      Bmat(2, pos) = -xiShape(1);

      pos = 0;
      AmatStress.setZero();
      prec sh, dsh;
      indexType cpos = pos + edgeNum * (this->intOrder + 1);
      triplets.clear();
      for (auto i = 0; i < this->intOrder + 1; ++i) {
        HierAMuS::LegendreShapes::getShape(sh, dsh, GP.getEta(ngp), i);
        indexType row = 0;
        triplets.emplace_back(row, cpos + i, sh);
      }
      pos += (this->intOrder + 1) * numberOfEdges;
      cpos = pos + this->intOrder * edgeNum;
      for (auto i = 0; i < this->intOrder; ++i) {
        HierAMuS::LegendreShapes::getShape(sh, dsh, GP.getEta(ngp), i);
        indexType row = 1;
        triplets.emplace_back(row, cpos + i, sh);
      }
      pos += (this->intOrder) * numberOfEdges;
      cpos = pos + this->intOrder * edgeNum;
      for (auto i = 0; i < this->intOrder; ++i) {
        HierAMuS::LegendreShapes::getShape(sh, dsh, GP.getEta(ngp), i);
        indexType row = 2;
        triplets.emplace_back(row, cpos + i, sh);
      }
      AmatSparseSig.setFromTriplets(triplets.begin(), triplets.end());

      AmatStress = AmatSparseSig;

      pos = 0;
      for (auto j = 2; j < WarpingShapeX.rows(); ++j) {
        AmatEpsC(0, pos) = WarpingShapeX(j);
        ++pos;
      }

      for (auto j = 1; j < dSurfShape.rows(); ++j) {
        AmatEpsC(1, pos) = dSurfShape(j);
        ++pos;
      }
      for (auto j = 2; j < WarpingShapeXY.rows(); ++j) {
        AmatEpsC(2, pos) = WarpingShapeXY(j);
        ++pos;
      }

      matData.strains.setZero();
      material->getMaterialData(pointers, matData, GP.getIntegrationPoint());

      Gmat += Bmat.transpose() * AmatStress * dA * GP.getWeight(ngp);
      Jmat += AmatStress.transpose() * matData.materialTangent * AmatStress *
              dA * GP.getWeight(ngp);
      Hmat +=
          AmatSparseSig.transpose() * AmatSparseSig * dA * GP.getWeight(ngp);

      Lmat += AmatStress.transpose() * AmatEpsC * dA * GP.getWeight(ngp);

      History.next();
    }
  }

}

void EL203_BeamInterface2D::setTangentResidualV6(
  PointerCollection& pointers,
  FiniteElement::beamInterfaceElement2D *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  Timer<HierAMuS::ms> testtimer;
  testtimer.start();
  auto GP = this->getIntegrationPoints(pointers, elem);
  indexType numberOfEdges = this->edgeMaterialMap.size();

  Types::VectorX<prec> WarpingShapeX, WarpingShapeXY, WarpingShapeY;
  Types::VectorX<prec> xiShape, dXiShape;
  Types::VectorX<prec> surfShape, dSurfShape;

  std::vector<DegreeOfFreedom *> surfaceDispDofs;
  elem->getDofsOnSolid(pointers, surfaceDispDofs, this->meshIdDisp,
                       this->intOrder);
  std::vector<DegreeOfFreedom *> vertexDispDofs;
  elem->getDofsOnVert(pointers, vertexDispDofs, this->meshIdDisp);
  std::vector<DegreeOfFreedom *> vertexRotDofs;
  elem->getDofsOnVert(pointers, vertexRotDofs, this->meshIdRot);

  Dofs.clear();
  Dofs.insert(Dofs.end(), surfaceDispDofs.begin(), surfaceDispDofs.end());
  Dofs.insert(Dofs.end(), vertexDispDofs.begin(), vertexDispDofs.end());
  Dofs.insert(Dofs.end(), vertexRotDofs.begin(), vertexRotDofs.end());

  //elem->getSolution(pointers, Dofs, solution);
  Types::VectorX<prec> solution = elem->getSolution(pointers, Dofs);
  

  indexType numDofs;
  numDofs = Dofs.size();
  stiffness.resize(numDofs, numDofs);
  residual.resize(numDofs);
  stiffness.setZero();
  residual.setZero();

  indexType numCols, numRows;

  numRows = 3;
  numCols = Dofs.size();
  Types::Matrix3X<prec> Bmat = Types::Matrix3X<prec>::Zero(numRows, numCols);

  numCols = (this->intOrder + 1) + this->intOrder * 2;
  numCols *= numberOfEdges;

  Types::Matrix3X<prec> AmatStress =
      Types::Matrix3X<prec>::Zero(numRows, numCols);
  Types::Matrix3X<prec> AmatTemp =
      Types::Matrix3X<prec>::Zero(numRows, numCols);


  Types::Matrix3X<prec> AmatEps = Types::Matrix3X<prec>::Zero(numRows, numCols);

  numCols = numberOfEdges * 3 * this->intOrder;

  Types::Matrix3X<prec> AmatEpsC =
      Types::Matrix3X<prec>::Zero(numRows, numCols - 2);

  Types::MatrixXX<prec> &Gmat = elem->get_element_data_field(pointers, m_G_id);
  Types::MatrixXX<prec> &Hmat = elem->get_element_data_field(pointers, m_H_id);
  Types::MatrixXX<prec> &Lmat = elem->get_element_data_field(pointers, m_L_id);
  Types::MatrixXX<prec> &Jmat = elem->get_element_data_field(pointers, m_J_id);

  Types::MatrixXX<prec> &resEpsC =
      elem->get_element_data_field(pointers, m_resEpsC_id);
  Types::MatrixXX<prec> &resEps =
      elem->get_element_data_field(pointers, m_resEps_id);
  Types::MatrixXX<prec> &resSig =
      elem->get_element_data_field(pointers, m_resSig_id);

  Types::MatrixXX<prec> &EpsC =
      elem->get_element_data_field(pointers, m_EpsC_id);
  Types::MatrixXX<prec> &Eps =
      elem->get_element_data_field(pointers, m_Eps_id);
  Types::MatrixXX<prec> &Sig =
      elem->get_element_data_field(pointers, m_Sig_id);
  

  Types::Vector3<prec> A1 = elem->getA1(pointers);
  Types::Vector3<prec> A2 = elem->getA2(pointers);

  std::vector<Eigen::Triplet<prec, indexType>> triplets;
  triplets.reserve((this->intOrder) * 3);
  Eigen::SparseMatrix<prec> AmatSparseSig(AmatStress.rows(), AmatStress.cols());

  testtimer.stop();
  std::cout << "Data allocation took: " << testtimer << std::endl;
  testtimer.start();

  auto History = elem->getHistoryDataIterator(pointers);

  Types::VectorX<prec> resEpsP(AmatEps.cols());
  Types::VectorX<prec> resSigP(AmatStress.cols());
  
  for (auto edgeNum = 0; edgeNum < numberOfEdges; ++edgeNum) {
    indexType gEdge = elem->getGlobalEdgeNumber(pointers, edgeNum);
    indexType mNum = this->edgeMaterialMap[gEdge];
    auto material =
        pointers.getMaterialFormulationList()->getMaterial(mNum);

    Materials::MaterialTransferData matData;
    matData.historyData = &History;

    for (auto ngp = 0; ngp < GP.getTotalGP(); ++ngp) {
      GP.setCurrNumber(ngp);
      prec dA =
          elem->getDA(pointers, edgeNum, GP.getXi(ngp), GP.getEta(ngp));
      Bmat.setZero();
      elem->getLocalWarpingShapesA1(pointers, WarpingShapeX,
                                    WarpingShapeXY, edgeNum, GP.getEta(ngp));
      elem->getXiShape(pointers, this->intOrder, GP.getXi(ngp), xiShape,
                       dXiShape);
      elem->getLocalSurfaceDispShapesSorted(pointers, surfShape,
                                            dSurfShape, edgeNum, GP.getEta(ngp),
                                            this->intOrder, this->meshIdDisp);
      prec z = elem->getZCoordinate(pointers, edgeNum, GP.getEta(ngp));

      indexType pos = 0;
      for (auto j = 0; j < surfaceDispDofs.size() / 3; ++j) {
        Bmat.block(0, pos, 1, 3) = A1.transpose() * dXiShape(0) * surfShape(j);
        Bmat.block(1, pos, 1, 3) = A2.transpose() * xiShape(0) * dSurfShape(j);
        Bmat.block(2, pos, 1, 3) = A1.transpose() * xiShape(0) * dSurfShape(j);
        Bmat.block(2, pos, 1, 3) += A2.transpose() * dXiShape(0) * surfShape(j);

        pos += 3;
      }
      Bmat.block(0, pos, 1, 3) = A1.transpose() * dXiShape(1);
      Bmat.block(2, pos, 1, 3) = A2.transpose() * dXiShape(1);
      pos += 3;
      Bmat(0, pos) = -z * dXiShape(1);
      Bmat(2, pos) = -xiShape(1);

      pos = 0;
      AmatStress.setZero();
      prec sh, dsh;
      indexType cpos = pos + edgeNum * (this->intOrder + 1);
      triplets.clear();
      for (auto i = 0; i < this->intOrder + 1; ++i) {
        HierAMuS::LegendreShapes::getShape(sh, dsh, GP.getEta(ngp), i);
        indexType row = 0;
        triplets.emplace_back(row, cpos + i, sh);
      }
      pos += (this->intOrder + 1) * numberOfEdges;
      cpos = pos + this->intOrder * edgeNum;
      for (auto i = 0; i < this->intOrder; ++i) {
        HierAMuS::LegendreShapes::getShape(sh, dsh, GP.getEta(ngp), i);
        indexType row = 1;
        triplets.emplace_back(row, cpos + i, sh);
      }
      pos += (this->intOrder) * numberOfEdges;
      cpos = pos + this->intOrder * edgeNum;
      for (auto i = 0; i < this->intOrder; ++i) {
        HierAMuS::LegendreShapes::getShape(sh, dsh, GP.getEta(ngp), i);
        indexType row = 2;
        triplets.emplace_back(row, cpos + i, sh);
      }
      AmatSparseSig.setFromTriplets(triplets.begin(), triplets.end());

      AmatStress = AmatSparseSig;

      pos = 0;
      for (auto j = 2; j < WarpingShapeX.rows(); ++j) {
        AmatEpsC(0, pos) = WarpingShapeX(j);
        ++pos;
      }

      for (auto j = 1; j < dSurfShape.rows(); ++j) {
        AmatEpsC(1, pos) = dSurfShape(j);
        ++pos;
      }
      for (auto j = 2; j < WarpingShapeXY.rows(); ++j) {
        AmatEpsC(2, pos) = WarpingShapeXY(j);
        ++pos;
      }

      matData.strains = Bmat * solution;
      material->getMaterialData(pointers, matData,GP.getIntegrationPoint());

      //Gmat += Bmat.transpose() * AmatStress * dA * GP.getWeight(ngp);
      Jmat += AmatStress.transpose() * matData.materialTangent * AmatStress *
              dA * GP.getWeight(ngp);
      //Hmat +=
      //    AmatSparseSig.transpose() * AmatSparseSig * dA * GP.getWeight(ngp);

      //Lmat += AmatStress.transpose() * AmatEpsC * dA * GP.getWeight(ngp);

      History.next();
    }
  }
  testtimer.stop();
  std::cout << "Integration duration: " << std::setprecision(16) << testtimer
            << " ms" << std::endl;

  // Eigen::EigenSolver<Types::MatrixXX<prec>> es(Jmat);
  // std::cout << "Eigenvalues of J are: " << es.eigenvalues() << std::endl;
  testtimer.start();
  // Types::MatrixXX<prec> HmatInv = Hmat.diagonal().asDiagonal().inverse();
  Types::VectorX<prec> HmatInv =
      Hmat.diagonal().asDiagonal().inverse().diagonal();

  Types::MatrixXX<prec> HMH =
      HmatInv.asDiagonal() * Jmat * HmatInv.asDiagonal();

  // Math::setEntriesToZeroEpsilon(HMH, prec(100));
  Types::MatrixXX<prec> LT = Lmat.transpose();
  Types::MatrixXX<prec> LHMHL = (LT * HMH * Lmat);
  // Math::setEntriesToZeroEpsilon(LHMHL, prec(100));
  Types::MatrixXX<prec> LHMHLLTinv = Math::AinvTimesB(LHMHL, LT);

  // stiffness = Gmat*HMH*Gmat.transpose() -
  // Gmat*HMH*Lmat*LHMHL*Lmat.transpose()*HMH*Gmat.transpose();
  Types::MatrixXX<prec> GHMH = Gmat * HMH;
  stiffness =
      GHMH * Gmat.transpose() - GHMH * Lmat * LHMHLLTinv * GHMH.transpose();

  testtimer.stop();

  std::cout << "First version duration: " << std::setprecision(16) << testtimer
            << " ms" << std::endl;



  residual = stiffness * solution;


}

void EL203_BeamInterface2D::setDegreesOfFreedomV7(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D &elem) {
  elem.setDofsOnSolid(pointers, this->meshIdDisp, this->intOrder);
  elem.setDofsOnVert(pointers, this->meshIdDisp);
  elem.setDofsOnVert(pointers, this->meshIdRot);
}

void EL203_BeamInterface2D::AdditionalOperationsV7(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D &elem) {
  elem.setBCOnAllNodesSolid(pointers, this->meshIdDisp, 2);

  // Restrict degrees of freedom on vertex
  elem.setBCOnVert(pointers, this->meshIdDisp, 2);
  elem.setBCOnVert(pointers, this->meshIdRot, 1);
  elem.setBCOnVert(pointers, this->meshIdRot, 2);

  elem.computeWarpingShapesNew(pointers, this->intOrder);
}

void EL203_BeamInterface2D::setTangentResidualV7(
  PointerCollection& pointers,
  FiniteElement::beamInterfaceElement2D *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  Timer<HierAMuS::ms> testtimer;
  testtimer.start();
  auto GP = this->getIntegrationPoints(pointers, elem);
  indexType numberOfEdges = this->edgeMaterialMap.size();

  Types::VectorX<prec> WarpingShapeX, WarpingShapeXY, WarpingShapeY,
      WarpingShapeYY;
  Types::VectorX<prec> solution, xiShape, dXiShape;
  Types::VectorX<prec> surfShape, dSurfShape;

  std::vector<DegreeOfFreedom *> surfaceDispDofs;
  elem->getDofsOnSolid(pointers, surfaceDispDofs, this->meshIdDisp,
                       this->intOrder);
  std::vector<DegreeOfFreedom *> vertexDispDofs;
  elem->getDofsOnVert(pointers, vertexDispDofs, this->meshIdDisp);
  std::vector<DegreeOfFreedom *> vertexRotDofs;
  elem->getDofsOnVert(pointers, vertexRotDofs, this->meshIdRot);

  Dofs.clear();
  Dofs.insert(Dofs.end(), surfaceDispDofs.begin(), surfaceDispDofs.end());
  Dofs.insert(Dofs.end(), vertexDispDofs.begin(), vertexDispDofs.end());
  Dofs.insert(Dofs.end(), vertexRotDofs.begin(), vertexRotDofs.end());

  elem->getSolution(pointers, Dofs, solution);

  indexType numDofs;
  numDofs = Dofs.size();
  stiffness.resize(numDofs, numDofs);
  residual.resize(numDofs);
  stiffness.setZero();
  residual.setZero();

  indexType numCols, numRows;

  Types::Matrix3X<prec> Bmat, AmatStress, AmatEps, AmatEpsC, AmatTemp;
  numRows = 3;
  numCols = Dofs.size();
  Bmat.resize(numRows, numCols);
  Bmat.setZero();

  numCols = (this->intOrder + 1) + this->intOrder * 2;
  numCols *= numberOfEdges;
  AmatStress.resize(numRows, numCols);
  AmatStress.setZero();
  AmatEps.resize(numRows, numCols);
  AmatEps.setZero();

  numCols = numberOfEdges * 3 * this->intOrder;
  indexType numCols2 = numberOfEdges * 2 * this->intOrder;
  AmatEpsC.resize(numRows, numCols2 - 1);
  AmatEpsC.setZero();

  AmatTemp = AmatStress;

  Types::MatrixXX<prec> Gmat, Hmat, Lmat, Jmat;
  Gmat.resize(Bmat.cols(), AmatStress.cols());
  Gmat.setZero();
  Hmat.resize(AmatStress.cols(), AmatEps.cols());
  Hmat.setZero();
  Lmat.resize(AmatStress.cols(), AmatEpsC.cols());
  Lmat.setZero();
  Jmat.resize(AmatEps.cols(), AmatEps.cols());
  Jmat.setZero();

  Types::Vector3<prec> A1;
  A1 = elem->getA1(pointers);
  Types::Vector3<prec> A2;
  A2 = elem->getA2(pointers);

  std::vector<Eigen::Triplet<prec, indexType>> triplets;
  triplets.reserve((this->intOrder) * 3);
  Eigen::SparseMatrix<prec> AmatSparseSig;
  AmatSparseSig.resize(AmatStress.rows(), AmatStress.cols());

  testtimer.stop();
  std::cout << "Data allocation took: " << testtimer << std::endl;
  testtimer.start();

  for (auto edgeNum = 0; edgeNum < numberOfEdges; ++edgeNum) {
    indexType gEdge = elem->getGlobalEdgeNumber(pointers, edgeNum);
    indexType mNum = this->edgeMaterialMap[gEdge];
    auto material =
        pointers.getMaterialFormulationList()->getMaterial(mNum);

    Materials::MaterialTransferData matData;

    for (auto ngp = 0; ngp < GP.getTotalGP(); ++ngp) {
      GP.setCurrNumber(ngp);
      prec dA =
          elem->getDA(pointers, edgeNum, GP.getXi(ngp), GP.getEta(ngp));
      Bmat.setZero();
      elem->getLocalWarpingShapesA1(pointers, WarpingShapeX,
                                    WarpingShapeXY, edgeNum, GP.getEta(ngp));
      elem->getLocalWarpingShapesA2(pointers, WarpingShapeY,
                                    WarpingShapeYY, edgeNum, GP.getEta(ngp));
      elem->getXiShape(pointers, this->intOrder, GP.getXi(ngp), xiShape,
                       dXiShape);
      elem->getLocalSurfaceDispShapesSorted(pointers, surfShape,
                                            dSurfShape, edgeNum, GP.getEta(ngp),
                                            this->intOrder, this->meshIdDisp);
      prec z = elem->getZCoordinate(pointers, edgeNum, GP.getEta(ngp));

      indexType pos = 0;
      for (auto j = 0; j < surfaceDispDofs.size() / 3; ++j) {
        Bmat.block(0, pos, 1, 3) = A1.transpose() * dXiShape(0) * surfShape(j);
        Bmat.block(1, pos, 1, 3) = A2.transpose() * xiShape(0) * dSurfShape(j);
        Bmat.block(2, pos, 1, 3) = A1.transpose() * xiShape(0) * dSurfShape(j);
        Bmat.block(2, pos, 1, 3) += A2.transpose() * dXiShape(0) * surfShape(j);
        // Bmat(0, pos) = dXiShape(0) * surfShape(j);
        // Bmat(1, pos + 1) = xiShape(0) * dSurfShape(j);
        // Bmat(2, pos) = xiShape(0) * dSurfShape(j);
        // Bmat(2, pos + 1) = dXiShape(0) * surfShape(j);

        pos += 3;
      }
      Bmat.block(0, pos, 1, 3) = A1.transpose() * dXiShape(1);
      Bmat.block(2, pos, 1, 3) = A2.transpose() * dXiShape(1);
      // Bmat(0, pos) = dXiShape(1);
      // Bmat(2, pos + 1) = dXiShape(1);
      pos += 3;
      Bmat(0, pos) = -z * dXiShape(1);
      Bmat(2, pos) = -xiShape(1);

      pos = 0;
      AmatStress.setZero();
      prec sh, dsh;
      indexType cpos = pos + edgeNum * (this->intOrder + 1);
      triplets.clear();
      for (auto i = 0; i < this->intOrder + 1; ++i) {
        HierAMuS::LegendreShapes::getShape(sh, dsh, GP.getEta(ngp), i);
        indexType row = 0;
        triplets.emplace_back(row, cpos + i, sh);
      }
      pos += (this->intOrder + 1) * numberOfEdges;
      cpos = pos + this->intOrder * edgeNum;
      for (auto i = 0; i < this->intOrder; ++i) {
        HierAMuS::LegendreShapes::getShape(sh, dsh, GP.getEta(ngp), i);
        indexType row = 1;
        triplets.emplace_back(row, cpos + i, sh);
      }
      pos += (this->intOrder) * numberOfEdges;
      cpos = pos + this->intOrder * edgeNum;
      for (auto i = 0; i < this->intOrder; ++i) {
        HierAMuS::LegendreShapes::getShape(sh, dsh, GP.getEta(ngp), i);
        indexType row = 2;
        triplets.emplace_back(row, cpos + i, sh);
      }
      AmatSparseSig.setFromTriplets(triplets.begin(), triplets.end());

      // AmatStress(0, pos + edgeNum) = prec(1);
      // pos += numberOfEdges;
      // AmatStress(1, pos + edgeNum) = prec(1);
      // pos += numberOfEdges;
      // for(auto i=0;i<numberOfEdges;++i){

      // AmatStress(1, pos) = dSurfShape(i);
      // pos += 1;
      // }
      // AmatStress(2, pos + edgeNum) = prec(1);
      AmatStress = AmatSparseSig;
      // AmatEps = AmatStress;

      pos = 0;
      for (auto j = 2; j < WarpingShapeX.rows(); ++j) {
        AmatEpsC(0, pos) = dXiShape(1) * WarpingShapeX(j);
        AmatEpsC(2, pos) = xiShape(1) * WarpingShapeXY(j);
        ++pos;
      }
      // AmatEpsC(1, pos + edgeNum) = prec(1);
      for (auto j = 1; j < WarpingShapeY.rows(); ++j) {
        AmatEpsC(1, pos) = xiShape(1) * WarpingShapeYY(j);
        AmatEpsC(2, pos) = dXiShape(1) * WarpingShapeY(j);
        ++pos;
      }
      for (auto j = 2; j < WarpingShapeXY.rows(); ++j) {
        // AmatEpsC(2, pos) = WarpingShapeXY(j);
        ++pos;
      }

      matData.strains = Bmat * solution;
      material->getMaterialData(pointers, matData,GP.getIntegrationPoint());

      // dA *= prec(2);

      Gmat += Bmat.transpose() * AmatStress * dA * GP.getWeight(ngp);
      //      Jmat += AmatEps.transpose() * matData.materialTangent * AmatEps *
      //      dA *
      //          GP.getWeight(ngp);
      Jmat += AmatStress.transpose() * matData.materialTangent * AmatStress *
              dA * GP.getWeight(ngp);
      // Hmat += AmatStress.transpose() * AmatEps * dA * GP.getWeight(ngp);
      Hmat +=
          AmatSparseSig.transpose() * AmatSparseSig * dA * GP.getWeight(ngp);

      Lmat += AmatStress.transpose() * AmatEpsC * dA * GP.getWeight(ngp);
    }
  }
  testtimer.stop();
  std::cout << "Integration duration: " << std::setprecision(16) << testtimer
            << " ms" << std::endl;

  // Eigen::EigenSolver<Types::MatrixXX<prec>> es(Jmat);
  // std::cout << "Eigenvalues of J are: " << es.eigenvalues() << std::endl;
  testtimer.start();
  // Types::MatrixXX<prec> HmatInv = Hmat.diagonal().asDiagonal().inverse();
  Types::VectorX<prec> HmatInv =
      Hmat.diagonal().asDiagonal().inverse().diagonal();

  Types::MatrixXX<prec> HMH =
      HmatInv.asDiagonal() * Jmat * HmatInv.asDiagonal();

  // Math::setEntriesToZeroEpsilon(HMH, prec(100));
  Types::MatrixXX<prec> LT = Lmat.transpose();
  Types::MatrixXX<prec> LHMHL = (LT * HMH * Lmat);
  // Math::setEntriesToZeroEpsilon(LHMHL, prec(100));
  Types::MatrixXX<prec> LHMHLLTinv = Math::AinvTimesB(LHMHL, LT);

  // stiffness = Gmat*HMH*Gmat.transpose() -
  // Gmat*HMH*Lmat*LHMHL*Lmat.transpose()*HMH*Gmat.transpose();
  Types::MatrixXX<prec> GHMH = Gmat * HMH;
  stiffness =
      GHMH * Gmat.transpose() - GHMH * Lmat * LHMHLLTinv * GHMH.transpose();

  testtimer.stop();

  std::cout << "First version duration: " << std::setprecision(16) << testtimer
            << " ms" << std::endl;

  // Types::MatrixXX<prec> JMatHmat, tempMat;
  // tempMat = Hmat.transpose();
  // JMatHmat = Math::AinvTimesB(Jmat, tempMat);

  // Types::MatrixXX<prec> Mmat = Hmat * JMatHmat;

  // Types::MatrixXX<prec> MinvG, MinvL;
  // tempMat = Gmat.transpose();
  // MinvG = Math::AinvTimesB(Mmat, tempMat);
  // MinvL = Math::AinvTimesB(Mmat, Lmat);

  // Types::MatrixXX<prec> GMinvL;
  // GMinvL = Gmat * MinvL;

  // Types::MatrixXX<prec> LMinvL;
  // LMinvL = Lmat.transpose() * MinvL;

  // Types::MatrixXX<prec> LMinvLInvGMinvL;
  // tempMat = GMinvL.transpose();
  // LMinvLInvGMinvL = Math::AinvTimesB(LMinvL, tempMat);

  // stiffness = Gmat * MinvG - GMinvL * LMinvLInvGMinvL;

  // auto Mmat = Hmat * Jmat.inverse() * Hmat.transpose();
  // auto Mmatinv = Mmat.inverse();
  // stiffness = Gmat * Mmatinv * Gmat.transpose() -
  //             Gmat * Mmatinv * Lmat *
  //                 (Lmat.transpose() * Mmatinv * Lmat).inverse() *
  //                 Lmat.transpose() * Mmatinv * Gmat.transpose();

  // Math::setEntriesToZeroEpsilon(stiffness, prec(100));

  residual = stiffness * solution;

  // indexType nn = Dofs.size();
  // nn=nn-6;
  // std::cout << stiffness.block(nn,nn,6,6) << std::endl;

  // Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

  // Eigen::EigenSolver<Types::MatrixXX<prec>> es(stiffness);
  // //std::cout << "Eigenvalues of J are:\n " << es.eigenvalues() << std::endl;
  // // std::cout << std::setprecision(3) << Hmat.format(CleanFmt)  <<
  // std::endl; elem->setUserConstant("eva", es.eigenvalues()(0).real());
  // elem->setUserConstant("evb", es.eigenvalues()(1).real());
  // elem->setUserConstant("evc", es.eigenvalues()(2).real());
}

auto EL203_BeamInterface2D::getNumberOfIntergrationPoints(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem)
    -> indexType {
  auto GP = this->getIntegrationPoints(pointers,elem);
  return GP.getTotalGP()*elem->getNumberOfEdges(pointers);
}

void EL203_BeamInterface2D::toParaviewAdaper(
    PointerCollection &pointers, FiniteElement::beamInterfaceElement2D &elem,
    vtkPlotInterface &paraviewAdapter, ParaviewSwitch control) {
#ifdef USE_VTK
  FiniteElement::beamInterfaceElement2D *element = &elem;

  if (this->plotVertexMap.empty()) {
    indexType nEdges = this->edgeMaterialMap.size();

    indexType maxVertNum = 0;
    for (auto i = 0; i < nEdges; ++i) {
      auto &edge =
          pointers.getGeometryData()->getEdgeData(element->getGlobalEdgeNumber(pointers, i));
      indexType nVerts = edge.getNumberOfVerts();
      for (auto j = 0; j < nVerts; ++j) {
        auto vert = edge.getVertex(j);
        indexType vertId = vert->getId();
        this->plotVertexMap[vertId] = j;
        if (vertId > maxVertNum)
          maxVertNum = vertId;
      }
    }
    ++maxVertNum;
    auto &vert = element->getVertex(pointers, 0);
    this->plotVertexMap[vert.getId()] = maxVertNum;
    if (vert.getId() > maxVertNum)
      maxVertNum = vert.getId();
    ++maxVertNum;
    for (auto i : this->plotVertexMap) {
      indexType id = i.first;
      this->plotVertexMapOpposite[id] = maxVertNum;
      ++maxVertNum;
    }
  }

  switch (control) {
  case ParaviewSwitch::Mesh: {
    indexType nEdges = this->edgeMaterialMap.size();

    indexType matNum = elem.getMaterial()->getNumber();
    std::vector<indexType> cellIds;
    for (auto i = 0; i < nEdges; ++i) {
      auto &edge =
          pointers.getGeometryData()->getEdgeData(element->getGlobalEdgeNumber(pointers, i));
      cellIds.resize(4);

      for (auto j = 0; j < 2; ++j) {
        auto vert = edge.getVertex(j);
        Types::Vector3<prec> coors = vert->getCoordinates();
        paraviewAdapter.addPoint(0, matNum, vert->getId(), coors);
        // std::cout << coors.transpose() << std::endl;
        coors += element->getA1(pointers) * element->getThickness(pointers);
        // std::cout << coors.transpose() << std::endl;
        paraviewAdapter.addPoint(
            0, matNum, this->plotVertexMapOpposite[vert->getId()], coors);
        cellIds[j * 3] = vert->getId();
        cellIds[1 + j % 2] = this->plotVertexMapOpposite[vert->getId()];
      }
      int celltype = VTK_QUAD;
      paraviewAdapter.addCell(0, matNum, elem.getId(), i + 1, cellIds, 4,
                              celltype);
    }
    auto &vert = element->getVertex(pointers, 0);
    Types::Vector3<prec> coors = vert.getCoordinates();
    paraviewAdapter.addPoint(0, matNum, vert.getId(), coors);
    cellIds.clear();
    cellIds.resize(1);
    cellIds[0] = vert.getId();
    paraviewAdapter.addCell(0, matNum, elem.getId(), nEdges + 1, cellIds, 1,
                            VTK_VERTEX);

  } break;
  case ParaviewSwitch::Solution: {
    std::vector<DegreeOfFreedom *> Dofs;
    Types::VectorX<prec> dispSolution, vertexDispSolution, vertexRotSolution;

    element->getDofsOnVert(pointers, Dofs, this->meshIdDisp);
    element->getSolution(pointers, Dofs, vertexDispSolution);
    Dofs.clear();
    element->getDofsOnVert(pointers, Dofs, this->meshIdRot);
    element->getSolution(pointers, Dofs, vertexRotSolution);

    indexType nEdges = this->edgeMaterialMap.size();

    std::vector<prec> sol, solr;
    sol.resize(3);
    solr.resize(3);

    Types::Vector3<prec> rsol, A1;
    A1 = element->getA1(pointers);

    indexType matNum = elem.getMaterial()->getNumber();
    for (auto i = 0; i < nEdges; ++i) {
      auto &edge =
          pointers.getGeometryData()->getEdgeData(element->getGlobalEdgeNumber(pointers, i));
      Dofs.clear();
      edge.getH1Dofs(Dofs, this->meshIdDisp, 1);
      element->getSolution(pointers, Dofs, dispSolution);
      prec eta = prec(-1);
      for (auto j = 0; j < 2; ++j) {
        auto vert = edge.getVertex(j);
        prec z = element->getZCoordinate(pointers, i, eta);
        eta += prec(2);
        rsol = vertexDispSolution - z * vertexRotSolution(0) * A1;

        for (auto k = 0; k < 3; ++k) {
          sol[k] = dispSolution(3 * j + k);
          solr[k] = rsol[k];
        }

        paraviewAdapter.setPointData(0, matNum, vert->getId(), sol, 3,
                                     paraviewNames::DisplacementName());
        paraviewAdapter.setPointData(
            0, matNum, this->plotVertexMapOpposite[vert->getId()], solr, 3,
            paraviewNames::DisplacementName());
      }
    }

    auto &vert = element->getVertex(pointers, 0);
    Dofs.clear();
    element->getDofsOnVert(pointers, Dofs, this->meshIdDisp);
    element->getSolution(pointers, Dofs, dispSolution);
    for (auto k = 0; k < 3; ++k) {
      sol[k] = dispSolution(k);
    }
    paraviewAdapter.setPointData(0, matNum, vert.getId(), sol, 3,
                                 paraviewNames::DisplacementName());

  } break;
  case HierAMuS::ParaviewSwitch::ProjectedValues: {

  } break;
  default:
    break;
  }

#endif
}

void EL203_BeamInterface2D::computeStressesAndStrains(
    FiniteElement::beamInterfaceElement2D &elem, indexType localEdgeNumber,
    prec xi, prec eta) {}

auto EL203_BeamInterface2D::getIntegrationPoints(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem)
    -> IntegrationPoints {
  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(this->intOrder * 2);
  return GP;
}




} // namespace HierAMuS::Elementformulations
