// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "beamInterfaceElement3D.h"

#include "MatrixTypes.h"
#include "datatypes.h"
#include "equations/DegreeOfFreedom.h"
#include "equations/Nodetypes.h"
#include "pointercollection/pointercollection.h"

#include "materials/MaterialformulationList.h"

#include "geometry/BeamInterface3D.h"
#include "geometry/GeometryData.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "shapefunctions/LegendreShapes.h"
#include "shapefunctions/LobattoShapes.h"
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include "control/HandlingStructs.h"
#include "control/OutputHandler.h"

#include "math/geometry.h"

#include <set>
#include <vtkCellType.h>

namespace HierAMuS::FiniteElement {
beamInterfaceElement3D::beamInterfaceElement3D() : m_vertMain(-1),m_vertX2(-1), m_vertX3(-1){}
beamInterfaceElement3D::~beamInterfaceElement3D() = default;

void beamInterfaceElement3D::setMaterialPerSubElement(
    std::vector<indexType> &materials) {
  m_materialPerFace = materials;
}

auto beamInterfaceElement3D::getMaterialFormulation(PointerCollection &pointers,
                                                    IntegrationPoint &ip)
    -> std::shared_ptr<HierAMuS::Materials::GenericMaterialFormulation> {
  return pointers.getMaterialFormulationList()->getMaterial(
      m_materialPerFace[ip.sectionNumber]);
}

void beamInterfaceElement3D::setH1Shapes(PointerCollection &pointers,
                                         indexType dispOrder,
                                         indexType meshIdDisp,
                                         indexType meshIdRot) {
  

  auto geometryData = pointers.getGeometryData();
  {
    auto face = geometryData->getFace(m_faces[0]);
    auto V1 = face->getVertex(pointers,0);
    minFaceVertNumber = V1->getId();
    maxFaceVertNumber = minFaceVertNumber;
  }
  for (auto i : m_faces) {
    auto face = geometryData->getFace(i);
    face->setH1Shapes(pointers, meshIdDisp, dispOrder, NodeTypes::displacement);
    for (auto j=0;j<face->getNumberOfVerts();++j)
    {
      auto V1 = face->getVertex(pointers, j);
      auto VId = V1->getId();
      if (minFaceVertNumber > VId)
        minFaceVertNumber = VId;
      if (maxFaceVertNumber < VId)
        maxFaceVertNumber = VId;
    }
  }
  auto &vert = geometryData->getVertex(m_beamVertex);
  vert.setNodeSet(pointers, meshIdDisp, 1, NodeTypes::displacement);
  vert.setNodeSet(pointers, meshIdRot, 1, NodeTypes::displacement);

  m_H1MeshId = meshIdDisp;
  m_MeshIdRot = meshIdRot;
  m_warpOrder = dispOrder;
}

auto beamInterfaceElement3D::getIntegrationPoints(ptrCol &pointers)
    -> IntegrationPoints {

  auto GP = pointers.getGeometryData()
                ->getFace(m_faces[0])
                ->getIntegrationPoints(pointers, this->id);
  GP.setType(IntegrationType::Scaled3D);
  GP.setNumberOfSections(m_faces.size());
  return GP;
}

void beamInterfaceElement3D::computeWarpingShapes(PointerCollection &pointers) {
  this->createNodeShapeMapping(pointers);
  this->computeGeometry(pointers);

  // displacement approximating shapes
  this->m_surfToBeamShapesU.resize(3, m_numberOfWarpingShapes * 3);
  this->m_surfToBeamShapesBeta.resize(3, m_numberOfWarpingShapes * 3);
  this->m_surfToBeamShapesU.setZero();
  this->m_surfToBeamShapesBeta.setZero();

  // real computation of the warping shape coefficients
  this->m_warpingCoefficients.resize(3, m_numberOfWarpingShapes);
  this->m_warpingCoefficients.setZero();

  if (m_warpingType == 3 || m_warpingType == 4) {
    m_warpingCoefficientsUy.resize(2, m_numberOfWarpingShapes);
    m_warpingCoefficientsUz.resize(2, m_numberOfWarpingShapes);
    m_warpingCoefficientsUy.setZero();
    m_warpingCoefficientsUz.setZero();
  }

  auto GP = this->getIntegrationPoints(pointers);
  GP.setOrder(m_warpOrder*2);

  Types::Matrix33<prec> eqSys;
  Types::Matrix33<prec> eqSysAverage;
  Types::Matrix22<prec> eqSys2;
  Types::Matrix22<prec> eqSys3;
  eqSysAverage.setZero();
  eqSys.setZero();
  eqSys2.setZero();
  eqSys3.setZero();

  prec A = prec(0);
  prec Ip = prec(0);
  prec I2 = prec(0);
  prec I3 = prec(0);
  prec I23 = prec(0);

  for (auto i : GP) {
    auto face = pointers.getGeometryData()->getFace(m_faces[i.sectionNumber]);
    auto Nodes = face->getH1Nodes(pointers, m_H1MeshId, m_warpOrder);
    auto jacobi = this->getJacobian(pointers, i);
    auto localCoord = this->getLocalCoordinate(pointers, i);
    IntegrationPoint faceip;
    faceip.xi = i.eta;
    faceip.eta = i.zeta;

    auto H1Shapes = face->getH1Shapes(pointers, m_warpOrder, faceip);

    auto dA = jacobi.block(1, 1, 2, 2).determinant() * i.weight;
    dA = abs(dA);
    eqSys(0, 0) += dA;
    eqSys(0, 1) += localCoord(1) * dA;
    eqSys(0, 2) += localCoord(2) * dA;
    eqSys(1, 1) += localCoord(1) * localCoord(1) * dA;
    eqSys(1, 2) += localCoord(1) * localCoord(2) * dA;
    eqSys(2, 2) += localCoord(2) * localCoord(2) * dA;

    A += dA;
    Ip += (localCoord(1) * localCoord(1) + localCoord(2) * localCoord(2)) * dA;
    I2 += (localCoord(2) * localCoord(2)) * dA;
    I3 += (localCoord(1) * localCoord(1)) * dA;
    I23 += (localCoord(1) * localCoord(2)) * dA;

    for (auto j = 0; j < Nodes.size(); ++j) {
      auto id = Nodes[j]->getId();
      indexType shapeNumber = this->m_nodeShapeMapping[id];
      this->m_warpingCoefficients(0, shapeNumber) -= H1Shapes.shapes(j) * dA;
      this->m_warpingCoefficients(1, shapeNumber) -=
          H1Shapes.shapes(j) * dA * localCoord(1);
      this->m_warpingCoefficients(2, shapeNumber) -=
          H1Shapes.shapes(j) * dA * localCoord(2);

      // New displacement approximating shapes
      this->m_surfToBeamShapesU(0, shapeNumber * 3) += H1Shapes.shapes(j) * dA;
      this->m_surfToBeamShapesU(1, shapeNumber * 3 + 1) +=
          H1Shapes.shapes(j) * dA;
      this->m_surfToBeamShapesU(2, shapeNumber * 3 + 2) +=
          H1Shapes.shapes(j) * dA;
      // New rotation approximating shapes
      this->m_surfToBeamShapesBeta(0, shapeNumber * 3 + 1) +=
          H1Shapes.shapes(j) * dA * (-localCoord(2));
      this->m_surfToBeamShapesBeta(0, shapeNumber * 3 + 2) +=
          H1Shapes.shapes(j) * dA * (localCoord(1));

      this->m_surfToBeamShapesBeta(1, shapeNumber * 3) +=
          H1Shapes.shapes(j) * dA * (localCoord(2));
      this->m_surfToBeamShapesBeta(2, shapeNumber * 3) +=
          H1Shapes.shapes(j) * dA * (-localCoord(1));
    }
  }

  eqSys(1, 0) = eqSys(0, 1);
  eqSys(2, 0) = eqSys(0, 2);
  eqSys(2, 1) = eqSys(1, 2);

  auto refVal = eqSys.norm();
  for (auto i = 0; i < 3; ++i) {
    for (auto j = 0; j < 3; ++j) {
      if (abs(eqSys(i, j)) < std::numeric_limits<prec>::epsilon() * refVal) {
        eqSys(i, j) = prec(0);
      }
    }
  }

  if (m_warpingType == 3 || m_warpingType == 4) {
    eqSys2(0, 0) = eqSys(0, 0);
    eqSys2(0, 1) = eqSys(0, 2);
    eqSys2(1, 0) = eqSys(2, 0);
    eqSys2(1, 1) = eqSys(2, 2);

    eqSys3(0, 0) = eqSys(0, 0);
    eqSys3(0, 1) = eqSys(0, 1);
    eqSys3(1, 0) = eqSys(1, 0);
    eqSys3(1, 1) = eqSys(1, 1);
    m_warpingCoefficientsUy.block(0, 0, 1, m_numberOfWarpingShapes) =
        m_warpingCoefficients.block(0, 0, 1, m_numberOfWarpingShapes);
    m_warpingCoefficientsUy.block(1, 0, 1, m_numberOfWarpingShapes) =
        m_warpingCoefficients.block(2, 0, 1, m_numberOfWarpingShapes);

    m_warpingCoefficientsUz.block(0, 0, 1, m_numberOfWarpingShapes) =
        m_warpingCoefficients.block(0, 0, 1, m_numberOfWarpingShapes);
    m_warpingCoefficientsUz.block(1, 0, 1, m_numberOfWarpingShapes) =
        m_warpingCoefficients.block(1, 0, 1, m_numberOfWarpingShapes);
    eqSys2 = eqSys2.inverse().eval();
    eqSys3 = eqSys3.inverse().eval();
    m_warpingCoefficientsUy = eqSys2 * m_warpingCoefficientsUy;
    m_warpingCoefficientsUz = eqSys3 * m_warpingCoefficientsUz;
  }

  eqSys = eqSys.inverse().eval();
  m_warpingCoefficients = eqSys * m_warpingCoefficients;

  // get verts with boundary conditions for warping shapes
  if (m_vertMain == -1 || m_vertX2 == -1 || m_vertX3 == -1)
    this->computeBCVerts(pointers);

  m_warpEqIds.resize(m_numberOfWarpingShapes * 3);
  for (auto i = 0; i < m_numberOfWarpingShapes * 3; ++i) {
    m_warpEqIds[i] = -1;
  }
  indexType totWarpDispFuncs = m_numberOfWarpingShapes * 3;
  { // First vertex, all Nodes are fixed
    auto &vt = pointers.getGeometryData()->getVertex(m_vertMain);
    auto NN = vt.getNodesOfSet(pointers, m_H1MeshId);
    indexType pos = m_nodeShapeMapping[NN[0]->getId()];
    m_warpEqIds[pos * 3] = totWarpDispFuncs - 1;
    m_warpEqIds[pos * 3 + 1] = totWarpDispFuncs - 2;
    m_warpEqIds[pos * 3 + 2] = totWarpDispFuncs - 3;
  }
  { // Second vertex, x and z dofs are fixed
    auto &vt = pointers.getGeometryData()->getVertex(m_vertX2);
    auto NN = vt.getNodesOfSet(pointers, m_H1MeshId);
    indexType pos = m_nodeShapeMapping[NN[0]->getId()];
    m_warpEqIds[pos * 3] = totWarpDispFuncs - 4;
    m_warpEqIds[pos * 3 + 2] = totWarpDispFuncs - 5;
  }
  { // Third vertex, x dof is fixed
    auto &vt = pointers.getGeometryData()->getVertex(m_vertX3);
    auto NN = vt.getNodesOfSet(pointers, m_H1MeshId);
    indexType pos = m_nodeShapeMapping[NN[0]->getId()];
    m_warpEqIds[pos * 3] = totWarpDispFuncs - 6;
  }
  indexType eqId = 0;
  for (auto &i:m_warpEqIds)
  {
    if (i==-1)
    {
      i = eqId;
      ++eqId;
    }
  }

  for (auto i : m_warpEqIds) {
    std::cout << i << " ";
  }
  std::cout << std::endl;

  m_surfToBeamShapesU = m_surfToBeamShapesU / A;
  eqSys.setZero();
  eqSys(0, 0) = Ip;
  eqSys(1, 1) = I2;
  eqSys(2, 2) = I3;
  eqSys(1, 2) = -I23;
  eqSys(2, 1) = -I23;

  eqSys = eqSys.inverse().eval();
  m_surfToBeamShapesBeta = (eqSys * m_surfToBeamShapesBeta).eval();
   
}

auto beamInterfaceElement3D::getWarpingShapes(
    PointerCollection &pointers, const Types::Matrix33<prec> &jacobi,
    IntegrationPoint &integrationPoint) -> warpingShapes3D {
  warpingShapes3D retVal;
  retVal.omega1.resize(m_numberOfWarpingShapes);
  retVal.omega2.resize(m_numberOfWarpingShapes);
  retVal.omega2.setZero();
  retVal.omega1Deriv.resize(3, m_numberOfWarpingShapes);
  retVal.omega2Deriv.resize(3, m_numberOfWarpingShapes);
  retVal.omega1Deriv.setZero();
  retVal.omega2Deriv.setZero();
  retVal.localShapes.resize(m_numberOfWarpingShapes);
  retVal.localShapes.setZero();
  retVal.localShapesDeriv.resize(3, m_numberOfWarpingShapes);
  retVal.localShapesDeriv.setZero();
  if (m_warpingType == 3 || m_warpingType == 4) {
    retVal.omega3.resize(this->m_numberOfWarpingShapes);
    retVal.omega3Deriv.resize(3, this->m_numberOfWarpingShapes);
    retVal.omega3.setZero();
    retVal.omega3Deriv.setZero();
  }

  auto face = pointers.getGeometryData()->getFace(
      m_faces[integrationPoint.sectionNumber]);
  auto Nodes = face->getH1Nodes(pointers, m_H1MeshId, m_warpOrder);
  IntegrationPoint faceIp;
  faceIp.xi = integrationPoint.eta;
  faceIp.eta = integrationPoint.zeta;
  auto localShapes = face->getH1Shapes(pointers, m_warpOrder, faceIp);

  Types::Matrix22<prec> jacobi2D = jacobi.block(1, 1, 2, 2);
  Types::Matrix22<prec> jacobi2Dinv = jacobi2D.inverse();
  localShapes.shapeDeriv = jacobi2Dinv.transpose() * localShapes.shapeDeriv;

  Types::Vector3<prec> localCoords =
      this->getLocalCoordinate(pointers, integrationPoint);

  retVal.omega1 =
      m_warpingCoefficients.block(0, 0, 1, m_numberOfWarpingShapes).transpose();
  retVal.omega1 += m_warpingCoefficients.block(1, 0, 1, m_numberOfWarpingShapes)
                       .transpose() *
                   localCoords(1);
  retVal.omega1 += m_warpingCoefficients.block(2, 0, 1, m_numberOfWarpingShapes)
                       .transpose() *
                   localCoords(2);

  retVal.omega1Deriv.block(1, 0, 1, m_numberOfWarpingShapes) =
      m_warpingCoefficients.block(1, 0, 1, m_numberOfWarpingShapes);
  retVal.omega1Deriv.block(2, 0, 1, m_numberOfWarpingShapes) =
      m_warpingCoefficients.block(2, 0, 1, m_numberOfWarpingShapes);

  if (m_warpingType == 3 || m_warpingType == 4) {
    retVal.omega2 =
        m_warpingCoefficientsUy.block(0, 0, 1, m_numberOfWarpingShapes)
            .transpose();
    retVal.omega2 +=
        m_warpingCoefficientsUy.block(1, 0, 1, m_numberOfWarpingShapes)
            .transpose() *
        localCoords(2);
    retVal.omega2Deriv.block(2, 0, 1, m_numberOfWarpingShapes) =
        m_warpingCoefficientsUy.block(1, 0, 1, m_numberOfWarpingShapes);

    retVal.omega3 =
        m_warpingCoefficientsUz.block(0, 0, 1, m_numberOfWarpingShapes)
            .transpose();
    retVal.omega3 +=
        m_warpingCoefficientsUz.block(1, 0, 1, m_numberOfWarpingShapes)
            .transpose() *
        localCoords(1);
    retVal.omega3Deriv.block(1, 0, 1, m_numberOfWarpingShapes) =
        m_warpingCoefficientsUz.block(1, 0, 1, m_numberOfWarpingShapes);
  }

  for (auto i = 0; i < Nodes.size(); ++i) {
    indexType id = Nodes[i]->getId();
    indexType shapeNumber = this->m_nodeShapeMapping[id];
    retVal.omega1(shapeNumber) += localShapes.shapes(i);
    retVal.omega2(shapeNumber) += localShapes.shapes(i);
    retVal.localShapes(shapeNumber) += localShapes.shapes(i);

    retVal.omega1Deriv(1, shapeNumber) += localShapes.shapeDeriv(0, i);
    retVal.omega1Deriv(2, shapeNumber) += localShapes.shapeDeriv(1, i);
    retVal.omega2Deriv(1, shapeNumber) += localShapes.shapeDeriv(0, i);
    retVal.omega2Deriv(2, shapeNumber) += localShapes.shapeDeriv(1, i);
    retVal.localShapesDeriv(1, shapeNumber) += localShapes.shapeDeriv(0, i);
    retVal.localShapesDeriv(2, shapeNumber) += localShapes.shapeDeriv(1, i);

    if (m_warpingType == 3 || m_warpingType == 4) {
      retVal.omega3(shapeNumber) += localShapes.shapes(i);
      retVal.omega3Deriv(1, shapeNumber) += localShapes.shapeDeriv(0, i);
      retVal.omega3Deriv(2, shapeNumber) += localShapes.shapeDeriv(1, i);
    }
  }

  // std::cout << retVal.omega1Deriv << std::endl;

  return retVal;
}

void beamInterfaceElement3D::setFaces(const std::vector<indexType> &facesIn) {
  m_faces = facesIn;
}

auto beamInterfaceElement3D::getH1Dofs(ptrCol &pointers, indexType meshIDdisp,
                                       indexType meshIDrot, indexType order)
    -> std::vector<DegreeOfFreedom *> {
  std::vector<DegreeOfFreedom *> Dofs;
  Dofs.resize(m_numberOfWarpingShapes * 3 + 6);
  for (auto i : m_faces) {
    auto face = pointers.getGeometryData()->getFace(i);
    auto nodes = face->getH1Nodes(pointers, meshIDdisp, order);
    for (auto j : nodes) {
      auto tempDof = j->getDegreesOfFreedom(pointers);
      indexType startpos = m_nodeShapeMapping[j->getId()] * 3;
      for (auto k : tempDof) {
        Dofs[startpos] = k;
        startpos++;
      }
    }
  }
  auto vert = pointers.getGeometryData()->getVertex(m_beamVertex);
  auto tn = vert.getNodesOfSet(pointers, meshIDdisp);
  auto tdof = tn[0]->getDegreesOfFreedom(pointers);
  indexType startpos = m_numberOfWarpingShapes * 3;
  for (auto j : tdof) {
    Dofs[startpos] = j;
    ++startpos;
  }
  tn = vert.getNodesOfSet(pointers, meshIDrot);
  tdof = tn[0]->getDegreesOfFreedom(pointers);
  for (auto j : tdof) {
    Dofs[startpos] = j;
    ++startpos;
  }
  return Dofs;
}

auto beamInterfaceElement3D::getRotationR0() -> Types::Matrix33<prec> {
  Types::Matrix33<prec> R0;
  R0.block(0, 0, 3, 1) = m_A1;
  R0.block(0, 1, 3, 1) = m_A2;
  R0.block(0, 2, 3, 1) = m_A3;
  return R0;
}

auto beamInterfaceElement3D::getNumberOfInnerDofs(indexType order)
    -> indexType {
  indexType op1 = order + 1;
  indexType e11 = op1 * op1;
  indexType e12 = op1 * order * 2;
  indexType t12 = op1 * order * 2;
  indexType t23 = order * order;

  indexType shapesPerNormalStress = e11 + e12;
  indexType shapesPerShearStrain = t12 + t23;
  indexType shapesPerFace = (shapesPerNormalStress + shapesPerShearStrain);

  indexType totShapes = shapesPerFace * this->m_faces.size();

  return totShapes;
}

auto beamInterfaceElement3D::getL2ShapesStressStrain(
    PointerCollection &pointers, IntegrationPoint &ip, indexType order,
    Types::MatrixXX<prec> &jacobi) -> localStressStrainInterface {

  localStressStrainInterface result;

  indexType shapesPerFace =
      this->getNumberOfInnerDofs(order) / this->m_faces.size();

  result.numShapes = shapesPerFace;
  result.pos = shapesPerFace * ip.sectionNumber;
  result.shapes.resize(6, shapesPerFace);
  result.shapes.setZero();

  Types::Matrix66<prec> transMat;
  transMat = this->getL2TransformationMatrix(pointers, ip);

  indexType pos = 0; // shapesPerFace*ip.sectionNumber;
  for (auto i = 0; i < 1; ++i) {
    auto LxShapes = LegendreShapes::getShape(ip.xi, i);
    for (auto j = 0; j <= order; ++j) {
      auto LeShapes = LegendreShapes::getShape(ip.eta, j);
      for (auto k = 0; k <= order; ++k) {
        auto LzShapes = LegendreShapes::getShape(ip.zeta, k);
        prec temp =
            LxShapes.shapeValue * LeShapes.shapeValue * LzShapes.shapeValue;
        if (i < 1) {
          result.shapes(0, pos) = temp;
          ++pos;
        }
        if (j < order) {
          result.shapes(1, pos) = temp;
          ++pos;
        }
        if (k < order) {
          result.shapes(2, pos) = temp;
          ++pos;
        }
        if (i < 1 && j < order) {
          result.shapes(3, pos) = temp;
          ++pos;
        }
        if (i < 1 && k < order) {
          result.shapes(4, pos) = temp;
          ++pos;
        }
        if (j < order && k < order) {
          result.shapes(5, pos) = temp;
          ++pos;
        }
      }
    }
  }

  // std::cout << "pos: " << pos << " max: " << shapesPerFace << std::endl;

  result.shapes = transMat * result.shapes;

  return result;
}

void beamInterfaceElement3D::setWarpingType(indexType typ) {
  m_warpingType = typ;
}

auto beamInterfaceElement3D::getNumberOfInnerWarpingDofs() -> indexType {

  indexType dofs;

  switch (m_warpingType) {
  case 1:
    dofs = this->getNumberOfInnerWarpingDofsV1();
    break;
  case 2:
    dofs = this->getNumberOfInnerWarpingDofsV2();
    break;
  case 3:
    dofs = this->getNumberOfInnerWarpingDofsV3();
    break;
  case 4:
    dofs = this->getNumberOfInnerWarpingDofsV4();
    break;
  }

  return dofs;
}

auto beamInterfaceElement3D::getWarpingMatrix(PointerCollection &pointers,
                                              IntegrationPoint &ip,
                                              indexType order,
                                              Types::MatrixXX<prec> &jacobi)
    -> Types::Matrix6X<prec> {
  switch (m_warpingType) {
  case 1:
    return this->getWarpingMatrixV1(pointers, ip, order, jacobi);
    break;
  case 2:
    return this->getWarpingMatrixV2(pointers, ip, order, jacobi);
    break;
  case 3:
    return this->getWarpingMatrixV3(pointers, ip, order, jacobi);
    break;
  case 4:
    return this->getWarpingMatrixV4(pointers, ip, order, jacobi);
    break;
  }

  return this->getWarpingMatrixV1(pointers, ip, order, jacobi);
}

void beamInterfaceElement3D::print(PointerCollection &pointers) {
  auto &Logger = pointers.getSPDLogger();

  Logger.info("Beam Interface Element 3D number");
  Logger.info("  H1 meshid:                 {:>12}", this->m_H1MeshId);
  Logger.info("  Warping order:             {:>12}", this->m_warpOrder);
  Logger.info("  Warping type:              {:>12}", this->m_warpingType);
  Logger.info("  Number of Warping shapes:  {:>12}", this->m_numberOfWarpingShapes);
  Logger.info("  Local Vector A1:           {}", this->m_A1.transpose());
  Logger.info("  Local Vector A2:           {}", this->m_A2.transpose());
  Logger.info("  Local Vector A3:           {}", this->m_A3.transpose());
  Logger.info("  Element thickness:         {:>12.6e}", this->m_length);


  Logger.debug("   shape function parameters");
  Logger.debug(this->m_warpingCoefficients.row(0));
  Logger.debug(this->m_warpingCoefficients.row(1));
  Logger.debug(this->m_warpingCoefficients.row(2));

}

auto beamInterfaceElement3D::getBMatrixWarpingLocal(
    PointerCollection &pointers, IntegrationPoint &ip, indexType meshId,
    indexType order, Types::MatrixXX<prec> &jacobi) -> Types::Matrix6X<prec> {
  auto face =
      pointers.getGeometryData()->getFace(this->m_faces[ip.sectionNumber]);
  auto Nodes = face->getH1Nodes(pointers, meshId, order);
  auto localCoord = this->getLocalCoordinate(pointers, ip);
  IntegrationPoint faceip;

  Types::Vector3<prec> normal = face->getFaceNormal(pointers);
  indexType beamShapeNumber = 1;
  prec beamFact = 1.0;
  if (normal.dot(m_A1) < 0.0) {
    beamShapeNumber = 0;
    beamFact = -1.0;
    std::cout << "Shape switch for beam side" << std::endl;
  }

  faceip.xi = ip.eta;
  faceip.eta = ip.zeta;
  auto localShapesT = face->getH1Shapes(pointers, m_warpOrder, faceip);

  auto sh = LobattoShapes::getShape(ip.xi, beamShapeNumber);

  localShapesT.shapeDeriv =
      jacobi.block<2, 2>(1, 1).inverse().transpose() * localShapesT.shapeDeriv;

  Types::Matrix6X<prec> B;
  B.resize(6, m_numberOfWarpingShapes * 3);
  B.setZero();
  sh.shapeDerivative *= (beamFact * prec(2) / m_length);

  Types::Matrix3X<prec> omh;
  omh.resize(3, m_surfToBeamShapesU.cols());
  omh = -m_surfToBeamShapesU;
  Types::Matrix3X<prec> omhxi2;
  omhxi2.resize(3, m_surfToBeamShapesU.cols());
  omhxi2.setZero();
  Types::Matrix3X<prec> omhxi3;
  omhxi3.resize(3, m_surfToBeamShapesU.cols());
  omhxi3.setZero();
  Types::Vector3<prec> rxi2;
  Types::Vector3<prec> rxi3;
  rxi2.setZero();
  rxi3.setZero();
  rxi2(1) = prec(1);
  rxi3(2) = prec(1);
  // Computing translational and rotational shape functions omh with derivatives
  localCoord(0) = prec(0);

 
  auto skewR = Math::Geometry::skewMatrix(localCoord);
  auto skewRxi2 = Math::Geometry::skewMatrix(rxi2);
  auto skewRxi3 = Math::Geometry::skewMatrix(rxi3);

  omh += skewR * m_surfToBeamShapesBeta;
  omhxi2 = skewRxi2 * m_surfToBeamShapesBeta;
  omhxi3 = skewRxi3 * m_surfToBeamShapesBeta;


  // Adding quadrilateral shape functions NS_hat
  for (indexType i = 0; i < Nodes.size(); ++i) {
    auto NN = Nodes[i];
    auto localId = m_nodeShapeMapping[NN->getId()];
    indexType pos = localId * 3;
    omh(0, pos + 0) += localShapesT.shapes(i);
    omh(1, pos + 1) += localShapesT.shapes(i);
    omh(2, pos + 2) += localShapesT.shapes(i);

    omhxi2(0, pos + 0) += localShapesT.shapeDeriv(0, i);
    omhxi2(1, pos + 1) += localShapesT.shapeDeriv(0, i);
    omhxi2(2, pos + 2) += localShapesT.shapeDeriv(0, i);

    omhxi3(0, pos + 0) += localShapesT.shapeDeriv(1, i);
    omhxi3(1, pos + 1) += localShapesT.shapeDeriv(1, i);
    omhxi3(2, pos + 2) += localShapesT.shapeDeriv(1, i);
  }

  // Apply zero translation and zero rotation to warp shape functions
  Types::Matrix6X<prec> condenMatrix;
  condenMatrix.resize(6, m_surfToBeamShapesU.cols());
  condenMatrix({0, 1, 2}, Eigen::all) = m_surfToBeamShapesU;
  condenMatrix({3, 4, 5}, Eigen::all) = m_surfToBeamShapesBeta;
  // sort condense matrix such that vb is at the end
  Types::Matrix6X<prec> condenMatrixSorted;
  condenMatrixSorted.resize(6, condenMatrix.cols());
  condenMatrixSorted(Eigen::all, m_warpEqIds) = condenMatrix;

  // extract block which needs to be inverted
  Types::Matrix66<prec> toInvert =
      condenMatrixSorted(Eigen::lastN(6), Eigen::lastN(6));
  Types::Matrix66<prec> Inverse = toInvert.inverse();

  // create modification matrix
  Types::Matrix6X<prec> modMatrix = Inverse * condenMatrixSorted;

  // modify omh
  //auto conden = [](Types::Matrix3X<prec> &toSort, std::vector<indexType> &ids,
  //                 Types::Matrix6X<prec> &mod) {
  //  using namespace Eigen;
  //  Types::Matrix3X<prec> Sorted;
  //  Sorted.resize(3, toSort.cols());
  //  Sorted(all, ids) = toSort;
  //  Types::Matrix3X<prec> temp2 = Sorted(all, lastN(6)) * mod;
  //  Types::Matrix3X<prec> temp = Sorted - temp2;
  //  toSort = temp(all, ids);
  //};

  //conden(omh, m_warpEqIds, modMatrix);
  //conden(omhxi2, m_warpEqIds, modMatrix);
  //conden(omhxi3, m_warpEqIds, modMatrix);

  // Adding rotational shape functions N_beta
  for (indexType i = 0; i < B.cols(); ++i) {
    B(0, i) = omh(0, i) * sh.shapeDerivative; // u_x,x
    B(1, i) = omhxi2(1, i) * sh.shapeValue;   // u_y,y
    B(2, i) = omhxi3(2, i) * sh.shapeValue;   // u_z,z
    B(3, i) = omh(1, i) * sh.shapeDerivative; // u_y,x
    B(3, i) = omhxi2(0, i) * sh.shapeValue;   // u_x,y
    B(4, i) = omh(2, i) * sh.shapeDerivative; // u_z,x
    B(4, i) = omhxi3(0, i) * sh.shapeValue;   // u_x,z
    B(5, i) = omhxi2(2, i) * sh.shapeValue;   // u_z,y
    B(5, i) = omhxi3(1, i) * sh.shapeValue;   // u_y,z
  }

  Types::Matrix6X<prec> B2;
  B2.resize(6, B.cols());
  B2(Eigen::all, m_warpEqIds) = B;
  Types::Matrix6X<prec> B3;
  B3.resize(6, B2.cols() - 6);
  B3.setZero();
  B3(Eigen::all, Eigen::all) =
      B2(Eigen::all, Eigen::seq(0, B3.cols() - 1));
  // B3(5, B3.cols() - 1) = (localCoord(1) + localCoord(2)) * sh.shapeValue;
  //B3(1, B3.cols() - 1) = sh.shapeValue;

  return B;
}

auto beamInterfaceElement3D::getBMatrixWarpingLocal2(
    PointerCollection &pointers, IntegrationPoint &ip, indexType meshId,
    indexType order, Types::MatrixXX<prec> &jacobi) -> Types::Matrix6X<prec> {

  Types::Matrix6X<prec> B(6,m_numberOfWarpingShapes*3+6);
  B.setZero();

  IntegrationPoint faceIp;
  faceIp.xi = ip.eta;
  faceIp.eta = ip.zeta;

  auto face = pointers.getGeometryData()->getFace(m_faces[ip.sectionNumber]);
  auto FaceShapes = face->getH1Shapes(pointers, order, faceIp);
  auto FaceNodes = face->getH1Nodes(pointers, m_H1MeshId, order);
  FaceShapes.shapeDeriv =
      jacobi.block<2, 2>(1, 1).inverse().transpose() * FaceShapes.shapeDeriv;

  for (indexType i = 0; i < FaceNodes.size();++i)
  {
    indexType pos = m_nodeShapeMapping[FaceNodes[i]->getId()] * 3;
    B.block<1,3>(1,pos)    = FaceShapes.shapeDeriv(0, i) * m_A2.transpose();
    B.block<1, 3>(2, pos)  = FaceShapes.shapeDeriv(1, i) * m_A3.transpose();
    B.block<1, 3>(3, pos)  = FaceShapes.shapeDeriv(0, i) * m_A1.transpose();
    B.block<1, 3>(4, pos)  = FaceShapes.shapeDeriv(1, i) * m_A1.transpose();
    B.block<1, 3>(5, pos)  = FaceShapes.shapeDeriv(1, i) * m_A2.transpose();
    B.block<1, 3>(5, pos) += FaceShapes.shapeDeriv(0, i) * m_A3.transpose();
  }
  //std::cout << "BMat local:\n" << B << std::endl;

  Types::Matrix3X<prec> omh(3, m_surfToBeamShapesU.cols());
  omh = m_surfToBeamShapesU;
  Types::Matrix3X<prec> omhxi2(3, m_surfToBeamShapesU.cols());
  omhxi2.setZero();
  Types::Matrix3X<prec> omhxi3(3, m_surfToBeamShapesU.cols());
  omhxi3.setZero();
  Types::Vector3<prec> rxi2;
  Types::Vector3<prec> rxi3;
  rxi2.setZero();
  rxi3.setZero();
  rxi2(1) = prec(1);
  rxi3(2) = prec(1);
  // Computing translational and rotational shape functions omh with derivatives
  auto localCoord = this->getLocalCoordinate(pointers, ip);

  auto skewR = Math::Geometry::skewMatrix(localCoord);
  auto skewRxi2 = Math::Geometry::skewMatrix(rxi2);
  auto skewRxi3 = Math::Geometry::skewMatrix(rxi3);

  omh += skewR * m_surfToBeamShapesBeta;
  omhxi2 = skewRxi2 * m_surfToBeamShapesBeta;
  omhxi3 = skewRxi3 * m_surfToBeamShapesBeta;

  auto sh1D = LobattoShapes::getShape(ip.xi, 1);
  sh1D.shapeDerivative *= (prec(2) / m_length);

  omh *= sh1D.shapeDerivative;
  omhxi2 *= sh1D.shapeValue;
  omhxi3 *= sh1D.shapeValue;

  Types::Matrix33<prec> R0 = this->getRotationR0().transpose();
  omh = R0 * omh;
  omhxi2 = R0 * omhxi2;
  omhxi3 = R0 * omhxi3;

  indexType rows = m_numberOfWarpingShapes * 3;
  B.block(0, 0, 1, rows) += omh.block(0, 0, 1, rows);
  B.block(1, 0, 1, rows) += omhxi2.block(1, 0, 1, rows);
  B.block(2, 0, 1, rows) += omhxi3.block(2, 0, 1, rows);

  B.block(3, 0, 1, rows) += omhxi2.block(0, 0, 1, rows);
  B.block(3, 0, 1, rows) += omh.block(1, 0, 1, rows);

  B.block(4, 0, 1, rows) += omhxi3.block(0, 0, 1, rows);
  B.block(4, 0, 1, rows) += omh.block(2, 0, 1, rows);

  B.block(5, 0, 1, rows) += omhxi3.block(1, 0, 1, rows);
  B.block(5, 0, 1, rows) += omhxi2.block(2, 0, 1, rows);


  indexType pos = rows;
  B.block(0, pos, 1, 3) = sh1D.shapeDerivative* m_A1.transpose();
  B.block(3, pos, 1, 3) = sh1D.shapeDerivative* m_A2.transpose();
  B.block(4, pos, 1, 3) = sh1D.shapeDerivative* m_A3.transpose();

  pos += 3;
  B.block(0, pos, 1, 3) =
      sh1D.shapeDerivative *
      (localCoord(2) * m_A2.transpose() - localCoord(1) * m_A3.transpose());

  B.block(3, pos, 1, 3) =
      -sh1D.shapeDerivative * localCoord(2) * m_A1.transpose();
  B.block(3, pos, 1, 3) += -sh1D.shapeValue * m_A3.transpose();

  B.block(4, pos, 1, 3) =
      sh1D.shapeDerivative * localCoord(1) * m_A1.transpose();
  B.block(4, pos, 1, 3) += sh1D.shapeValue * m_A2.transpose();

  return B;
}

auto beamInterfaceElement3D::getWarpDofConstraint(PointerCollection &pointers)
    -> std::set<indexType> {
  // Computing columns to remove from B
  auto &v1 = pointers.getGeometryData()->getVertex(m_vertMain);
  auto &v2 = pointers.getGeometryData()->getVertex(m_vertX2);
  auto &v3 = pointers.getGeometryData()->getVertex(m_vertX3);

  auto N1 = v1.getNodesOfSet(pointers, m_H1MeshId);
  auto N2 = v2.getNodesOfSet(pointers, m_H1MeshId);
  auto N3 = v3.getNodesOfSet(pointers, m_H1MeshId);
  std::set<indexType> removeCols;
  // Main vertex columns
  indexType pos1 = m_nodeShapeMapping[N1[0]->getId()] * 3;
  removeCols.insert(pos1);
  removeCols.insert(pos1 + 1);
  removeCols.insert(pos1 + 2);
  // xi2 vertex columns
  indexType pos2 = m_nodeShapeMapping[N2[0]->getId()] * 3;
  removeCols.insert(pos2);
  removeCols.insert(pos2 + 2);
  // xi2 vertex columns
  indexType pos3 = m_nodeShapeMapping[N3[0]->getId()] * 3;
  removeCols.insert(pos3);
  return removeCols;
}

void beamInterfaceElement3D::setWarpDofConstraintVertices(
    std::vector<indexType> &vertnumbers)
{
  m_vertMain = vertnumbers[0];
  m_vertX2 = vertnumbers[1];
  m_vertX3 = vertnumbers[2];
}

void beamInterfaceElement3D::geometryToParaview(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh)
{
  indexType numPoints = 8;
  std::vector<indexType> points(numPoints);
  indexType cc=1;
  for (auto ff : m_faces)
  {
    auto face = pointers.getGeometryData()->getFace(ff);
    indexType nVerts = face->getNumberOfVerts();
    for (indexType nn=0;nn<nVerts;++nn)
    {
      auto V1 = face->getVertex(pointers, nn);
      V1->geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
      Geometry::Vertex V2;
      V2.setCoordinates(V1->getCoordinates() + m_A1 * m_length);
      V2.setId(V1->getId() + maxFaceVertNumber);
      points[nn] = V1->getId();
      points[nn + 4] = V2.getId();
      V1->geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
      V2.geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
    }
    paraviewAdapter.addCell(mainMesh, subMesh, this->id, cc, points, numPoints,
                            VTK_HEXAHEDRON);
    cc++;
  }
}

void beamInterfaceElement3D::computeWeightsParaview(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh) {}

void beamInterfaceElement3D::H1SolutionToParaview(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType meshId, indexType order,
    std::string name)
{
  auto &BeamVert = pointers.getGeometryData()->getVertex(m_beamVertex);
  auto beamDispNodes = BeamVert.getNodesOfSet(pointers, m_H1MeshId);
  auto beamRotNodes = BeamVert.getNodesOfSet(pointers, m_MeshIdRot);
  auto beamDispDofs = beamDispNodes[0]->getDegreesOfFreedom(pointers);
  auto beamRotDofs = beamRotNodes[0]->getDegreesOfFreedom(pointers);
  Types::Vector3<prec> beamDispSol =
      pointers.getSolutionState()->getSolution(beamDispDofs);
  Types::Vector3<prec> beamRotSol =
      pointers.getSolutionState()->getSolution(beamRotDofs);

  //auto beamCoor = BeamVert.getCoordinates();
  auto R0 = this->getRotationR0().transpose();

  auto allDofs = this->getH1Dofs(pointers, m_H1MeshId, m_MeshIdRot, m_warpOrder);
  allDofs.erase(allDofs.end() - 6, allDofs.end());

  Types::VectorX<prec> allSol =
      pointers.getSolutionState()->getSolution(allDofs);


  for (auto nf : m_faces)
  {
    auto face = pointers.getGeometryData()->getFace(nf);
    auto fDofs = face->getH1Dofs(pointers, m_H1MeshId, 1);
    Types::VectorX<prec> fSol = pointers.getSolutionState()->getSolution(fDofs);
    for (auto nV=0;nV<face->getNumberOfVerts();++nV)
    {
      auto V1 = face->getVertex(pointers, nV);
      Types::Vector3<prec> cCoor = V1->getCoordinates() - m_projectedCoordinate;
      cCoor = R0 * cCoor;
      indexType idA = V1->getId();
      indexType idB = idA + maxFaceVertNumber;
      auto VNodes = V1->getNodesOfSet(pointers, m_H1MeshId);
      Types::Vector3<prec> VSol =
          pointers.getSolutionState()->getSolution(*VNodes[0]);
      std::vector<prec> solVec(3);
      for (auto i=0;i<3;++i)
      {
        solVec[i] = VSol(i);
      }
      paraviewAdapter.setPointData(mainMesh, subMesh, idA, solVec, 3,
                                   paraviewNames::DisplacementName());

      auto skewR = Math::Geometry::skewMatrix(cCoor);
      Types::Vector3<prec> wtSol =
          (-m_surfToBeamShapesU + skewR * m_surfToBeamShapesBeta) * allSol;
      //indexType pos = m_nodeShapeMapping[VNodes[0]->getId()] * 3;
      VSol += beamDispSol - skewR * beamRotSol + wtSol;
      for (auto i = 0; i < 3; ++i) {
        solVec[i] = VSol(i);
      }
      paraviewAdapter.setPointData(mainMesh, subMesh, idB, solVec, 3,
                                   paraviewNames::DisplacementName());
    }
  }
}

void beamInterfaceElement3D::projectDataToParaviewVertices(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {}

auto beamInterfaceElement3D::getBeamShapesSpecial(PointerCollection &pointers,
                                                  IntegrationPoint &ip,
                                                  indexType order,
                                                  Types::MatrixXX<prec> &jacobi)
    -> beamShapesSpecial {
  return {};
}

void beamInterfaceElement3D::setVerts(std::vector<indexType> &vertsIn) {
  m_beamVertex = vertsIn[0];
}

auto beamInterfaceElement3D::getElement(PointerCollection &pointers)
    -> Geometry::BeamInterface3D * {
  auto elem = pointers.getGeometryData()->getSpecial(m_geoNumber);
  return dynamic_cast<Geometry::BeamInterface3D *>(elem);
}

void beamInterfaceElement3D::createNodeShapeMapping(
    PointerCollection &pointers) {
  indexType ns = 0;
  for (auto faceId : m_faces) {
    auto face = pointers.getGeometryData()->getFace(faceId);
    for (auto i = 0; i < face->getNumberOfVerts(); ++i) {
      auto Vert = face->getVertex(pointers, i);
      std::vector<GenericNodes *> Nodes;
      Vert->getNodes(pointers, Nodes, m_H1MeshId);
      for (auto &node : Nodes) {
        if (m_nodeShapeMapping.find(node->getId()) ==
            m_nodeShapeMapping.end()) {
          m_nodeShapeMapping[node->getId()] = ns;
          ++ns;
        }
      }
    }
  }
  for (auto faceId : m_faces) {
    auto face = pointers.getGeometryData()->getFace(faceId);
    auto Nodes = face->getH1Nodes(pointers, m_H1MeshId, m_warpOrder);
    // auto Vert=  face->getVertex(pointers, 0);
    for (auto &node : Nodes) {
      if (m_nodeShapeMapping.find(node->getId()) == m_nodeShapeMapping.end()) {
        m_nodeShapeMapping[node->getId()] = ns;
        ++ns;
      }
    }
  }
  m_numberOfWarpingShapes = ns;
}

void beamInterfaceElement3D::computeGeometry(PointerCollection &pointers) {
  auto face = pointers.getGeometryData()->getFace(m_faces[0]);
  Types::Vector3<prec> x1;
  Types::Vector3<prec> x2;
  Types::Vector3<prec> x3;
  x1 = face->getVertex(pointers, 0)->getCoordinates();
  x2 = face->getVertex(pointers, 1)->getCoordinates();
  x3 = face->getVertex(pointers, 2)->getCoordinates();

  // compute the local basis system on the surface.
  m_A2 = x2 - x1;
  m_A3 = x3 - x1;
  m_A1 = m_A2.cross(m_A3).normalized();
  m_A2 = m_A3.cross(m_A1).normalized();
  m_A3 = m_A1.cross(m_A2).normalized();
  m_A1 = {prec(1), 0, 0};
  m_A2 = {0, prec(1), 0};
  m_A3 = {0, 0, prec(1)};

  // compute the projected coordinate.
  Types::Vector3<prec> beam_coordinate =
      pointers.getGeometryData()->getVertex(m_beamVertex).getCoordinates();
  Types::Vector3<prec> temp = beam_coordinate - x1;
  m_length = temp.dot(m_A1);

  m_shapeXiFace = 0;
  m_shapeXiBeam = 1;

  if (m_length < 0) {
    m_A1 = -m_A1;
    // m_A2 = -m_A2;
    m_A3 = (m_A1.cross(m_A2)).normalized();
    m_length = -m_length;
    m_shapeXiFace = 1;
    m_shapeXiBeam = 0;
  }

  m_projectedCoordinate = beam_coordinate - m_length * m_A1;
}

void beamInterfaceElement3D::computeFunctionParameters(
    PointerCollection &pointers) {
  // Loop over all faces
  //for (auto i : this->m_faces) {
  //  auto *face = pointers.getGeometryData()->getFace(i);
  //  auto Nodes = face->getH1Nodes(pointers, m_H1MeshId, m_warpOrder);
  //  auto GP = face->getIntegrationPoints(pointers, -1);
  //  //indexType totShapes = m_numberOfWarpingShapes;
//
  //  GP.setOrder(m_warpOrder);
  //  for (auto i : GP) {
  //    auto shapes = face->getH1Shapes(pointers, m_warpOrder, i);
  //    for (auto j = 0; j < Nodes.size(); ++j) {
  //      indexType nId = Nodes[j]->getId();
  //    }
  //  }
  //}
}

auto beamInterfaceElement3D::getL2TransformationMatrix(
    PointerCollection &pointers, IntegrationPoint &ip)
    -> Types::Matrix66<prec> {
  IntegrationPoint centerIp;
  centerIp.xi = prec(0);
  centerIp.eta = prec(0);
  centerIp.zeta = prec(0);
  centerIp.sectionNumber = ip.sectionNumber;

  Types::Matrix33<prec> J0 = this->getJacobian(pointers, centerIp);

  Types::Matrix66<prec> trans;
  trans(0, 0) = J0(0, 0) * J0(0, 0);
  trans(0, 1) = J0(1, 0) * J0(1, 0);
  trans(0, 2) = J0(2, 0) * J0(2, 0);
  trans(0, 3) = J0(0, 0) * J0(1, 0);
  trans(0, 4) = J0(0, 0) * J0(2, 0);
  trans(0, 5) = J0(1, 0) * J0(2, 0);
  trans(1, 0) = J0(0, 1) * J0(0, 1);
  trans(1, 1) = J0(1, 1) * J0(1, 1);
  trans(1, 2) = J0(2, 1) * J0(2, 1);
  trans(1, 3) = J0(0, 1) * J0(1, 1);
  trans(1, 4) = J0(0, 1) * J0(2, 1);
  trans(1, 5) = J0(1, 1) * J0(2, 1);
  trans(2, 0) = J0(0, 2) * J0(0, 2);
  trans(2, 1) = J0(1, 2) * J0(1, 2);
  trans(2, 2) = J0(2, 2) * J0(2, 2);
  trans(2, 3) = J0(0, 2) * J0(1, 2);
  trans(2, 4) = J0(0, 2) * J0(2, 2);
  trans(2, 5) = J0(1, 2) * J0(2, 2);
  trans(3, 0) = 2 * J0(0, 0) * J0(0, 1);
  trans(3, 1) = 2 * J0(1, 0) * J0(1, 1);
  trans(3, 2) = 2 * J0(2, 0) * J0(2, 1);
  trans(3, 3) = J0(0, 0) * J0(1, 1) + J0(0, 1) * J0(1, 0);
  trans(3, 4) = J0(0, 0) * J0(2, 1) + J0(0, 1) * J0(2, 0);
  trans(3, 5) = J0(1, 0) * J0(2, 1) + J0(1, 1) * J0(2, 0);
  trans(4, 0) = 2 * J0(0, 0) * J0(0, 2);
  trans(4, 1) = 2 * J0(1, 0) * J0(1, 2);
  trans(4, 2) = 2 * J0(2, 0) * J0(2, 2);
  trans(4, 3) = J0(0, 0) * J0(1, 2) + J0(0, 2) * J0(1, 0);
  trans(4, 4) = J0(0, 0) * J0(2, 2) + J0(0, 2) * J0(2, 0);
  trans(4, 5) = J0(1, 0) * J0(2, 2) + J0(1, 2) * J0(2, 0);
  trans(5, 0) = 2 * J0(0, 1) * J0(0, 2);
  trans(5, 1) = 2 * J0(1, 1) * J0(1, 2);
  trans(5, 2) = 2 * J0(2, 1) * J0(2, 2);
  trans(5, 3) = J0(0, 1) * J0(1, 2) + J0(0, 2) * J0(1, 1);
  trans(5, 4) = J0(0, 1) * J0(2, 2) + J0(0, 2) * J0(2, 1);
  trans(5, 5) = J0(1, 1) * J0(2, 2) + J0(1, 2) * J0(2, 1);

  return trans;
}

auto beamInterfaceElement3D::getLocalCoordinate(
    PointerCollection &pointers, IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  auto face = pointers.getGeometryData()->getFace(
      m_faces[integrationPoint.sectionNumber]);
  IntegrationPoint faceInt;
  faceInt.xi = integrationPoint.eta;
  faceInt.eta = integrationPoint.zeta;
  Types::Vector3<prec> coor = face->getCoordinates(pointers, faceInt);
  coor -= m_projectedCoordinate;
  coor = this->getRotationR0().transpose() * coor;
  return coor;
}

auto beamInterfaceElement3D::getJacobian(PointerCollection &pointers,
                                         IntegrationPoint &integrationPoint)
    -> Types::MatrixXX<prec> {
  auto face = pointers.getGeometryData()->getFace(
      m_faces[integrationPoint.sectionNumber]);
  IntegrationPoint faceInt;
  faceInt.xi = integrationPoint.eta;
  faceInt.eta = integrationPoint.zeta;
  Types::MatrixXX<prec> jacobi(3, 3);
  Types::Vector3<prec> normal = face->getFaceNormal(pointers);
  jacobi.block(0, 0, 3, 1) = normal * m_length * prec(0.5);
  jacobi.block(0, 1, 3, 1) = face->getTangent_G1(pointers, faceInt);
  jacobi.block(0, 2, 3, 1) = face->getTangent_G2(pointers, faceInt);

  jacobi = getRotationR0().transpose() * jacobi;

  if (jacobi.determinant() < 0.0) {
    throw std::runtime_error("Error: Negative Jacobian");
  }

  return jacobi;
}

auto beamInterfaceElement3D::getH1Shapes(PointerCollection &pointers,
                                         indexType dispOrder,
                                         indexType meshIdDisp,
                                         const Types::Matrix33<prec> &jacobi,
                                         IntegrationPoint &integrationPoint)
    -> Geometry::H1Shapes {
  auto face = pointers.getGeometryData()->getFace(
      m_faces[integrationPoint.sectionNumber]);

  indexType totshapes = m_numberOfWarpingShapes + 2;
  IntegrationPoint faceIntegrationPoint;
  faceIntegrationPoint.xi = integrationPoint.eta;
  faceIntegrationPoint.eta = integrationPoint.zeta;

  Types::Vector3<prec> normal = face->getFaceNormal(pointers);
  indexType faceShapeNumber = 0;
  indexType beamShapeNumber = 1;
  prec beamFact = 1.0;
  if (normal.dot(m_A1) < 0.0) {
    faceShapeNumber = 1;
    beamShapeNumber = 0;
    beamFact = -1.0;
  }

  auto localShapes =
      face->getH1Shapes(pointers, dispOrder, faceIntegrationPoint);
  auto faceNodes = face->getH1Nodes(pointers, meshIdDisp, dispOrder);

  auto shx = LobattoShapes::getShape(integrationPoint.xi, faceShapeNumber);

  Geometry::H1Shapes shapes;
  shapes.shapes.resize(totshapes);
  shapes.shapeDeriv.resize(3, totshapes);
  shapes.shapes.setZero();
  shapes.shapeDeriv.setZero();

  for (auto i = 0; i < localShapes.shapes.rows(); ++i) {
    indexType shindex = m_nodeShapeMapping[faceNodes[i]->getId()];
    shapes.shapes(shindex) = localShapes.shapes(i) * shx.shapeValue;
    shapes.shapeDeriv(0, shindex) = localShapes.shapes(i) * shx.shapeDerivative;
    shapes.shapeDeriv(1, shindex) =
        localShapes.shapeDeriv(0, i) * shx.shapeValue;
    shapes.shapeDeriv(2, shindex) =
        localShapes.shapeDeriv(1, i) * shx.shapeValue;
  }
  shapes.shapeDeriv = jacobi.inverse().transpose() * shapes.shapeDeriv;

  shx = LobattoShapes::getShape(integrationPoint.xi, beamShapeNumber);

  shapes.shapes(totshapes - 2) = shx.shapeValue;
  shapes.shapes(totshapes - 1) = shx.shapeValue;
  shapes.shapeDeriv(0, totshapes - 2) =
      shx.shapeDerivative * prec(2) / m_length * beamFact;
  shapes.shapeDeriv(0, totshapes - 1) =
      shx.shapeDerivative * prec(2) / m_length * beamFact;

  return shapes;
}

auto beamInterfaceElement3D::getNumberOfInnerWarpingDofsV1() -> indexType {
  indexType dofs;

  dofs = m_numberOfWarpingShapes * 4;
  dofs -= 3; // Ex
  dofs -= 3; // Exy
  dofs -= 2; // Ey
  dofs -= 2; // Ez

  return dofs;
}

auto beamInterfaceElement3D::getWarpingMatrixV1(PointerCollection &pointers,
                                                IntegrationPoint &ip,
                                                indexType order,
                                                Types::MatrixXX<prec> &jacobi)
    -> Types::Matrix6X<prec> {
  indexType numShapes = this->getNumberOfInnerWarpingDofsV1();
  Types::Matrix6X<prec> result(6, numShapes);
  result.setZero();

  auto shapes = this->getWarpingShapes(pointers, jacobi, ip);
  indexType pos = 0;
  // std::cout << "omega1: \n" << shapes.omega1 << std::endl;
  // std::cout << "omega1Deriv: \n" << shapes.omega1Deriv << std::endl;
  for (auto i = 3; i < m_numberOfWarpingShapes; ++i) {
    result(0, pos) = shapes.omega1(i); // Exx
    ++pos;
    result(3, pos) = shapes.omega1Deriv(1, i); // tauxy
    result(4, pos) = shapes.omega1Deriv(2, i); // tauxz
    ++pos;
  }
  auto face = pointers.getGeometryData()->getFace(this->m_faces[0]);
  auto V1 = face->getVertex(pointers, 0);
  auto V2 = face->getVertex(pointers, 1);
  auto coor = V2->getCoordinates() - V1->getCoordinates();
  auto dotpV = coor.dot(m_A2);

  // Case parallel to a2
  indexType start = 1;
  indexType deriva = 1;
  indexType derivb = 2;
  indexType rem = 2;
  // Case parallel to a3
  if (abs(dotpV) < std::numeric_limits<prec>::epsilon() * prec(10000)) {
    start = 2;
    deriva = 2;
    derivb = 1;
    rem = 1;
  }
  for (auto i = 1; i < m_numberOfWarpingShapes; ++i) {
    if (i != rem) {
      result(start, pos) = shapes.omega2Deriv(deriva, i);
      result(5, pos) = shapes.omega2Deriv(derivb, i);
      ++pos;
    }
  }
  // Case parallel to a3
  if (start == 2) {
    start = 1; // Eyy
    deriva = 1;
    derivb = 2;
    rem = 2;
  } else { // Case parallel to a2
    start = 2;
    deriva = 2;
    derivb = 1;
    rem = 2;
  }
  for (auto i = 1; i < m_numberOfWarpingShapes; ++i) {
    if (i != rem) {
      result(start, pos) = shapes.omega2Deriv(deriva, i);
      result(5, pos) = shapes.omega2Deriv(derivb, i);
      ++pos;
    }
  }
  // auto hh = this->getLocalCoordinate(pointers, ip);
  // result(5, pos) = hh(1) * hh(2);
  // std::cout << pos << " " << numShapes << std::endl;
  return result;
}

auto beamInterfaceElement3D::getNumberOfInnerWarpingDofsV2() -> indexType {
  indexType dofs;

  dofs = m_numberOfWarpingShapes * 4;
  dofs -= 3; // Ex
  dofs -= 3; // Exy
  dofs -= 2; // Ey
  dofs -= 2; // Ez
  dofs -= 0; // Eyz

  return dofs;
}

auto beamInterfaceElement3D::getWarpingMatrixV2(PointerCollection &pointers,
                                                IntegrationPoint &ip,
                                                indexType order,
                                                Types::MatrixXX<prec> &jacobi)
    -> Types::Matrix6X<prec> {
  indexType numShapes = this->getNumberOfInnerWarpingDofsV2();
  Types::Matrix6X<prec> result(6, numShapes);
  result.setZero();

  auto shapes = this->getWarpingShapes(pointers, jacobi, ip);
  indexType pos = 0;
  // std::cout << "omega1: \n" << shapes.omega1 << std::endl;
  // std::cout << "omega1Deriv: \n" << shapes.omega1Deriv << std::endl;
  for (auto i = 3; i < m_numberOfWarpingShapes; ++i) {
    result(0, pos) = shapes.omega1(i); // Exx
    ++pos;
    result(3, pos) = shapes.omega1Deriv(1, i); // tauxy
    result(4, pos) = shapes.omega1Deriv(2, i); // tauxz
    ++pos;
  }
  auto face = pointers.getGeometryData()->getFace(this->m_faces[0]);
  auto V1 = face->getVertex(pointers, 0);
  auto V2 = face->getVertex(pointers, 1);
  auto coor = V2->getCoordinates() - V1->getCoordinates();
  auto dotpV = coor.dot(m_A2);

  // Case parallel to a2
  indexType start = 1;
  indexType deriva = 1;
  indexType rem = 2;
  // Case parallel to a3
  if (abs(dotpV) < std::numeric_limits<prec>::epsilon() * prec(10000)) {
    start = 2;
    deriva = 2;
    rem = 1;
  }
  for (auto i = 1; i < m_numberOfWarpingShapes; ++i) {
    if (i != rem) {
      result(start, pos) = shapes.omega2Deriv(deriva, i);
      // result(5, pos) = shapes.omega2Deriv(derivb, i);
      ++pos;
    }
  }
  // Case parallel to a3
  if (start == 2) {
    start = 1; // Eyy
    deriva = 1;
    rem = 2;
  } else { // Case parallel to a2
    start = 2;
    deriva = 2;
    rem = 1;
  }
  for (auto i = 1; i < m_numberOfWarpingShapes; ++i) {
    if (i != rem) {
      result(start, pos) = shapes.omega2Deriv(deriva, i);
      // result(5, pos) = shapes.omega2Deriv(derivb, i);
      ++pos;
    }
  }

  // for (auto i = 1; i < m_numberOfWarpingShapes; ++i) {
  //   //if (i != rem) {
  //     //result(start, pos) = shapes.omega2Deriv(deriva, i);
  //   result(5, pos) =
  //       shapes.omega2Deriv(deriva, i);
  //   ++pos;
  //   result(5, pos) =
  //       shapes.omega2Deriv(derivb, i);
  //   ++pos;
  //   //}
  // }
  //  auto hh = this->getLocalCoordinate(pointers, ip);
  //  result(5, pos) = hh(1) * hh(2);
  // std::cout << pos << " " << numShapes << std::endl;
  return result;
}

auto beamInterfaceElement3D::getNumberOfInnerWarpingDofsV3() -> indexType {
  indexType dofs;

  dofs = m_numberOfWarpingShapes * 4;
  dofs -= 3; // Ex
  dofs -= 3; // Exy
  dofs -= 2; // Ey
  dofs -= 2; // Ez
  dofs -= 0; // Eyz

  return dofs;
}

auto beamInterfaceElement3D::getWarpingMatrixV3(PointerCollection &pointers,
                                                IntegrationPoint &ip,
                                                indexType order,
                                                Types::MatrixXX<prec> &jacobi)
    -> Types::Matrix6X<prec> {
  indexType numShapes = this->getNumberOfInnerWarpingDofsV3();
  Types::Matrix6X<prec> result(6, numShapes);
  result.setZero();

  auto shapes = this->getWarpingShapes(pointers, jacobi, ip);
  indexType pos = 0;
  // std::cout << "omega1: \n" << shapes.omega1 << std::endl;
  // std::cout << "omega1Deriv: \n" << shapes.omega1Deriv << std::endl;
  for (auto i = 3; i < m_numberOfWarpingShapes; ++i) {
    result(0, pos) = shapes.omega1(i); // Exx
    ++pos;
    result(3, pos) = shapes.omega1Deriv(1, i); // tauxy
    //++pos;
    result(4, pos) = shapes.omega1Deriv(2, i); // tauxz
    ++pos;
  }
  auto face = pointers.getGeometryData()->getFace(this->m_faces[0]);
  auto V1 = face->getVertex(pointers, 0);
  auto V2 = face->getVertex(pointers, 1);
  auto coor = V2->getCoordinates() - V1->getCoordinates();
  auto dotpV = coor.dot(m_A2);

  // Case parallel to a2
  indexType rem = 2;
  // Case parallel to a3
  if (abs(dotpV) < std::numeric_limits<prec>::epsilon() * prec(10000)) {
    rem = 1;
  }
  for (auto i = 1; i < m_numberOfWarpingShapes; ++i) {
    if (i != rem) {
      result(1, pos) = shapes.omega2Deriv(1, i);
      //++pos;
      result(5, pos) = shapes.omega2Deriv(2, i);
      ++pos;
    }
  }
  // std::cout << shapes.omega2Deriv << "\n" << std::endl;
  //  Case parallel to a3
  if (rem == 2) {
    rem = 1;
  } else { // Case parallel to a2
    rem = 2;
  }
  for (auto i = 1; i < m_numberOfWarpingShapes; ++i) {
    if (i != rem) {
      result(2, pos) = shapes.omega3Deriv(2, i);
      //++pos;
      result(5, pos) = shapes.omega3Deriv(1, i);
      ++pos;
    }
  }
  // std::cout << " test " << std::endl;
  // std::cout << numShapes << " " << pos << "\n" << std::endl;

  // for (auto i = 1; i < m_numberOfWarpingShapes; ++i) {
  //   //if (i != rem) {
  //     //result(start, pos) = shapes.omega2Deriv(deriva, i);
  //   result(5, pos) =
  //       shapes.omega2Deriv(deriva, i);
  //   ++pos;
  //   result(5, pos) =
  //       shapes.omega2Deriv(derivb, i);
  //   ++pos;
  //   //}
  // }
  //  auto hh = this->getLocalCoordinate(pointers, ip);
  //  result(5, pos) = hh(1) * hh(2);
  return result;
}

auto beamInterfaceElement3D::getNumberOfInnerWarpingDofsV4() -> indexType {
  indexType dofs;

  dofs = m_numberOfWarpingShapes * 6;
  dofs -= 3; // Ex
  dofs -= 2; // Exy
  dofs -= 2; // Exy
  dofs -= 1; // Ey
  dofs -= 1; // Ez
  dofs -= 1; // Eyz

  return dofs;
}

auto beamInterfaceElement3D::getWarpingMatrixV4(PointerCollection &pointers,
                                                IntegrationPoint &ip,
                                                indexType order,
                                                Types::MatrixXX<prec> &jacobi)
    -> Types::Matrix6X<prec> {
  indexType numShapes = this->getNumberOfInnerWarpingDofsV4();
  Types::Matrix6X<prec> result(6, numShapes);
  result.setZero();

  auto shapes = this->getWarpingShapes(pointers, jacobi, ip);
  indexType pos = 0;

  // Eps_x
  for (auto i = 3; i < m_numberOfWarpingShapes; ++i) {
    result(0, pos) = shapes.omega1(i); // Exx
    ++pos;
  }

  for (auto i = 1; i < m_numberOfWarpingShapes; ++i) {
    // Eps_y
    result(1, pos) = shapes.localShapesDeriv(1, i);
    ++pos;
    // Eps_z
    result(2, pos) = shapes.localShapesDeriv(2, i);
    ++pos;
    // Eps_yz
    result(5, pos) =
        shapes.localShapesDeriv(2, i) + shapes.localShapesDeriv(1, i);
    ++pos;
  }
  auto face = pointers.getGeometryData()->getFace(this->m_faces[0]);
  auto V1 = face->getVertex(pointers, 0);
  auto V2 = face->getVertex(pointers, 1);
  auto coor = V2->getCoordinates() - V1->getCoordinates();
  auto dotpV = coor.dot(m_A2);

  // Case parallel to a2
  indexType start = 2;
  indexType rem = 1;
  // Case parallel to a3
  if (abs(dotpV) < std::numeric_limits<prec>::epsilon() * prec(10000)) {
    start = 1;
    rem = 2;
  }
  // Eps_xy
  for (auto i = start; i < m_numberOfWarpingShapes; ++i) {
    if (i != rem) {
      result(3, pos) = shapes.omega2Deriv(1, i);
      ++pos;
    }
  }
  // Case parallel to a3
  if (start == 2) {
    start = 1;
    rem = 2;
  } else { // Case parallel to a2
    start = 2;
    rem = 1;
  }
  for (auto i = start; i < m_numberOfWarpingShapes; ++i) {
    if (i != rem) {
      result(4, pos) = shapes.omega3Deriv(2, i);
      ++pos;
    }
  }

  return result;
}
void beamInterfaceElement3D::computeBCVerts(PointerCollection &pointers) {

  auto R = this->getRotationR0();
  std::set<indexType> vertNums;
  for (auto i : m_faces) {
    auto face = pointers.getGeometryData()->getFace(i);
    std::vector<indexType> lvertnums;
    face->getVerts(lvertnums);
    for (auto nn : lvertnums) {
      vertNums.insert(nn);
    }
  }

  struct vnumCoor {
    indexType num;
    Types::Vector3<prec> coor;

    void setFromOther(vnumCoor &other, Types::Vector3<prec> &refcoor) {
      prec normOther = (other.coor - refcoor).norm();
      prec normThis = (this->coor - refcoor).norm();
      if (normOther < normThis) {
        this->num = other.num;
        this->coor = other.coor;
      }
    }
  };

  prec ymin = std::numeric_limits<prec>::max();
  prec ymax = std::numeric_limits<prec>::min();
  prec zmin = ymin;
  prec zmax = ymax;

  for (auto vnum : vertNums) {
    auto &vert = pointers.getGeometryData()->getVertex(vnum);
    Types::Vector3<prec> coor =
        R.transpose() * (vert.getCoordinates() - m_projectedCoordinate);
    if (coor(1) > ymax)
      ymax = coor(1);
    if (coor(1) < ymin)
      ymin = coor(1);
    if (coor(2) > zmax)
      zmax = coor(2);
    if (coor(2) < zmin)
      zmin = coor(2);
  }

  Types::Vector3<prec> c1, c2, c3;
  c1.setZero();
  c2.setZero();
  c3.setZero();

  c1(1) = ymin;
  c1(2) = zmin;
  c2 = c1;
  c3 = c1;
  c2(1) = ymax;
  c3(2) = zmax;

  vnumCoor v1, v2, v3;
  {
    indexType vv = *vertNums.begin();
    v1.num = vv;
    Types::Vector3<prec> tt =
        pointers.getGeometryData()->getVertex(vv).getCoordinates();
    v1.coor = R.transpose() * (tt - m_projectedCoordinate);
  }
  v2 = v1;
  v3 = v1;

  for (auto vnum : vertNums) {
    auto &vert = pointers.getGeometryData()->getVertex(vnum);
    vnumCoor currV;
    currV.num = vnum;
    currV.coor =
        R.transpose() * (vert.getCoordinates() - m_projectedCoordinate);

    v1.setFromOther(currV, c1);
    v2.setFromOther(currV, c2);
    v3.setFromOther(currV, c3);
  }

  std::cout << "num: " << v1.num << " coor: " << v1.coor.transpose()
            << std::endl;
  std::cout << "num: " << v2.num << " coor: " << v2.coor.transpose()
            << std::endl;
  std::cout << "num: " << v3.num << " coor: " << v3.coor.transpose()
            << std::endl;

  m_vertMain = v1.num;
  m_vertX2 = v2.num;
  m_vertX3 = v3.num;
}
} // namespace HierAMuS::FiniteElement
