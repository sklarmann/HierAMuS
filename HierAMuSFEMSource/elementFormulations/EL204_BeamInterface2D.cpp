// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include <iostream>


#include <elementFormulations/EL204_BeamInterface2D.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include "control/ParameterList.h"
#include <pointercollection/pointercollection.h>

#include <finiteElements/GenericFiniteElement.h>
#include <finiteElements/beamInterfaceElement2D.h>


#include "shapefunctions/IntegrationsPoints/dataClasses/GaussPoints.h"

#include <elementFormulations/GenericElementFormulation.h>
#include <solver/GenericSolutionState.h>

#include <types/MatrixTypes.h>

namespace HierAMuS::Elementformulations {

EL204_BeamInterface2D::EL204_BeamInterface2D(PointerCollection *ptrCol)
    : GenericElementFormulationInterface(ptrCol) {}

EL204_BeamInterface2D::~EL204_BeamInterface2D() = default;

void EL204_BeamInterface2D::readData(PointerCollection &pointers,
                                     ParameterList &list) {
  this->meshIdDisp = list.getIndexVal("meshiddisp");
  this->meshIdWarp = list.getIndexVal("meshidwarp");
  this->meshIdRot = list.getIndexVal("meshidrot");
  this->planeStrain = list.getIndexVal("planestrain");
  this->mode = list.getIndexVal("mode");

  this->intOrder = list.getIndexVal("intorder");
  this->emod = list.getPrecVal("emodul");
  this->nu = list.getPrecVal("nu");
  this->thick = list.getPrecVal("thick");


  auto &Log = pointers.getSPDLogger();

  Log.info("\n{:-<100}\n"
                "*   Element 204, specified Options\n"
                "    Mesh id for displacement nodes:         {:>12}\n"
                "    Shape function order for displacements: {:>12}\n"
                "{:-<100}\n",
                "",
                this->meshIdDisp,
                this->intOrder,
                "");


  this->messageUnprocessed(pointers, list, "EL204_BeamInterface2D");
}

void EL204_BeamInterface2D::setDegreesOfFreedom(
    PointerCollection &pointers, FiniteElement::beamInterfaceElement2D &elem) {

  auto *ee = dynamic_cast<FiniteElement::beamInterfaceElement2D *>(&elem);
  // ee->setShapes(*this->ptrCol, this->meshIdDisp, this->meshIdWarp,
  // this->meshIdRot, this->intOrder);

  ee->setDofsOnSolid(pointers, this->meshIdDisp, this->intOrder);
  ee->setDofsOnSolid(pointers, this->meshIdWarp, this->intOrder);
  ee->setDofsOnVert(pointers, this->meshIdDisp);
  ee->setDofsOnVert(pointers, this->meshIdRot);

  ee->setMeshIdDisp(this->meshIdDisp);
  ee->setMeshIdRot(this->meshIdRot);
  ee->setMeshIdWarp(this->meshIdWarp);
}


auto EL204_BeamInterface2D::getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
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


void EL204_BeamInterface2D::AdditionalOperations(
  PointerCollection& pointers, FiniteElement::beamInterfaceElement2D &elem) {
  auto *ee = dynamic_cast<FiniteElement::beamInterfaceElement2D *>(&elem);

  // Geometry Set up
  ee->computeGeometry(pointers);

  // Degree of Freedom set up
  ee->setNodeMapDisp(pointers, this->meshIdDisp);
  ee->setNodeMapWarp(pointers, this->meshIdWarp);

  // Restrict unused nodes on Solid
  ee->setBCOnAllNodesSolid(pointers, this->meshIdDisp, 2);
  ee->setBCOnAllNodesSolid(pointers, this->meshIdWarp, 2);
  ee->setBCOnAllNodesSolid(pointers, this->meshIdWarp, 1);
  // ee->setBCOnAllNodesSolid(*this->ptrCol, this->meshIdWarp, 0);
  indexType numSolidNodes =
      ee->getNumberOfSolidNodes(pointers, this->meshIdDisp);
  ee->setBCOnNodeSolid(pointers, this->meshIdWarp, 0, 0);
  ee->setBCOnNodeSolid(pointers, this->meshIdWarp, 0, 1);
  ee->setBCOnNodeSolid(pointers, this->meshIdWarp, 0, 2);
  ee->setBCOnNodeSolid(pointers, this->meshIdWarp, numSolidNodes - 1, 0);

  // Restrict degrees of freedom on vertex
  ee->setBCOnVert(pointers, this->meshIdDisp, 2);
  ee->setBCOnVert(pointers, this->meshIdRot, 1);
  ee->setBCOnVert(pointers, this->meshIdRot, 2);

  ee->computeShapesLocalWarping(pointers, this->meshIdWarp,
                                this->intOrder);
  ee->computeSurfaceDispShapes(pointers, this->meshIdDisp, this->intOrder);
}

void EL204_BeamInterface2D::setTangentResidual(
    PointerCollection &pointers, FiniteElement::beamInterfaceElement2D &elem,
  Types::MatrixXX<prec> &stiffness, Types::VectorX<prec> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  auto *ee = dynamic_cast<FiniteElement::beamInterfaceElement2D *>(&elem);
  indexType numEdges = ee->getNumberOfEdges(pointers);

  Types::Matrix33<prec> cmat;
  cmat.setZero();

  if (this->planeStrain == 1) {
    prec fac =
        this->emod / ((prec)1 + this->nu) / ((prec)1 - (prec)2 * this->nu);
    cmat(0, 0) = ((prec)1 - this->nu);
    cmat(0, 1) = (this->nu);
    cmat(1, 1) = ((prec)1 - this->nu);
    cmat(1, 0) = (this->nu);
    cmat *= fac;
    cmat(2, 2) = this->emod / ((prec)1 + this->nu) / (prec)2;

  } else {
    prec fac = this->emod / ((prec)1 - this->nu * this->nu);
    cmat(0, 0) = ((prec)1);
    cmat(0, 1) = (this->nu);
    cmat(1, 1) = ((prec)1);
    cmat(1, 0) = (this->nu);
    cmat *= fac;
    cmat(2, 2) = this->emod / ((prec)1 + this->nu) / (prec)2;
    cmat *= this->thick;
  }

  IntegrationPoints GP = IntegrationPointsManagement::getIntegrationsPoints(-1);
  GP.setTypeOrder(IntegrationType::Gauss2D, this->intOrder + 1);

  Types::VectorX<prec> WarpingShapeX, WarpingShapeXY, WarpingShapeY,
      WarpingShapeYY, Nu, Nbeta, xiShape, dXiShape, surfShape, dSurfShape;
  Types::VectorX<prec> tempConst, tempLinear;

  ee->getSurfaceDispShapes(Nu, Nbeta);

  // std::cout << Nu << std::endl;
  // std::cout << Nbeta << std::endl;

  std::vector<DegreeOfFreedom *> surfaceDispDofs;
  ee->getDofsOnSolid(pointers, surfaceDispDofs, this->meshIdDisp,
                     this->intOrder);
  std::vector<DegreeOfFreedom *> surfaceWarpDofs;
  ee->getDofsOnSolid(pointers, surfaceWarpDofs, this->meshIdWarp,
                     this->intOrder);
  std::vector<DegreeOfFreedom *> vertexDispDofs;
  ee->getDofsOnVert(pointers, vertexDispDofs, this->meshIdDisp);
  std::vector<DegreeOfFreedom *> vertexRotDofs;
  ee->getDofsOnVert(pointers, vertexRotDofs, this->meshIdRot);

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
  //prec Asum = 0, Asum2 = 0;

  for (auto i = 0; i < numEdges; ++i) {
    prec dA = ee->getJacXi() * ee->getJacEta(pointers, i);
    //Asum += dA;
    for (auto ngp = 0; ngp < GP.getTotalGP(); ++ngp) {
      Bmat.setZero();
      ee->getWarpingShapes(pointers, WarpingShapeX, WarpingShapeY,
                           WarpingShapeXY, WarpingShapeYY, i, GP.getEta(ngp),
                           this->intOrder, this->meshIdWarp);

      ee->getXiShape(pointers, this->intOrder, GP.getXi(ngp), xiShape,
                     dXiShape);
      ee->getLocalSurfaceDispShapesSorted(pointers, surfShape, dSurfShape,
                                          i, GP.getEta(ngp), this->intOrder,
                                          this->meshIdDisp);
      auto z = ee->getZCoordinate(pointers, i, GP.getEta(ngp));

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

      // std::cout << GP.getXi(ngp) << " " << GP.getEta(ngp) <<
      // GP.getWeight(ngp) << std::endl; std::cout << Bmat << std::endl;
      //Asum2 += dA * GP.getWeight(ngp);
      stiffness += Bmat.transpose() * cmat * Bmat * dA * GP.getWeight(ngp);
    }
  }
  
}

auto EL204_BeamInterface2D::getNumberOfIntergrationPoints(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem)
    -> indexType {
  IntegrationPoints GP = IntegrationPointsManagement::getIntegrationsPoints(-1);
  GP.setTypeOrder(IntegrationType::Gauss2D, this->intOrder + 1);
  auto *ee = dynamic_cast<FiniteElement::beamInterfaceElement2D *>(elem);
  indexType numEdges = ee->getNumberOfEdges(pointers);
  return GP.getTotalGP()*numEdges;
}

void EL204_BeamInterface2D::toParaviewAdaper(
    PointerCollection &pointers, FiniteElement::beamInterfaceElement2D &elem,
    vtkPlotInterface &paraviewAdapter, ParaviewSwitch control)
{
  pointers.getSPDLogger().warn("Plot functionality for EL204_BeamInterface2D not implemented!");

}

void EL204_BeamInterface2D::setTangentResidualMixed(
    FiniteElement::GenericFiniteElement *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {}


} // namespace HierAMuS

