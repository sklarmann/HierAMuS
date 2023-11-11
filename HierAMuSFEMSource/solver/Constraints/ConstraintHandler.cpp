// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "ConstraintHandler.h"
#include "GeneralLink.h"
#include "MatrixTypes.h"
#include "control/BinaryWrite.h"
#include "datatypes.h"
#include "math/EigenSparseOperations.h"

// Geometry
#include "geometry/GeometryData.h"
#include "geometry/Edges/LinearEdgeRuntime.h"
#include "geometry/Faces/FacesRuntime.h"


#include "control/BinaryWrite.h"

#include "Timer.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "geometry/VertexData.h"
#include "geometry/Edges/EdgesData.h"
#include "geometry/Faces/FacesData.h"

// Equations
#include "EquationHandler.h"

namespace HierAMuS {
ConstraintHandler::ConstraintHandler() {}

ConstraintHandler::ConstraintHandler(const ConstraintHandler &other) {
  this->m_A = other.m_A;
  this->m_ATranspose = other.m_ATranspose;
  this->m_dB = other.m_dB;
  this->m_NewtonInc = other.m_NewtonInc;

  this->m_constraints.reserve(other.m_constraints.size());
  for (auto i : other.m_constraints) {
    this->m_constraints.push_back(i->getCopy());
  }
}

void ConstraintHandler::initialize(PointerCollection &pointers) {
  auto eqHandler = pointers.getEquationHandler();
  auto activeIds = eqHandler->getNumberOfActiveEquations();
  auto inActiveIds = eqHandler->getNumberOfInActiveEquations();
  auto totalIds = eqHandler->getNumberOfTotalEquations();

  m_A.resize(inActiveIds, activeIds);
  m_dB.resize(inActiveIds);
  m_NewtonInc.resize(totalIds);

  Types::VectorX<indexType> entriesPerCol;
  entriesPerCol.resize(activeIds);
  entriesPerCol.setZero();

  for (auto Constr : m_constraints) {
    auto DofR = Constr->getDofRelation(pointers);
    for (auto j : DofR.m_slaveDofs) {

      m_NewtonInc.coeffRef(j->getId()) = prec(1);
      m_dB.coeffRef(j->getEqId()) = prec(1);
    }
    for (auto j : DofR.m_masterDofs) {
      auto slaveId = j->getEqId();
      entriesPerCol(slaveId) += DofR.m_masterDofs.size();
    }
  }

  setEigenSparseVectorZero(m_dB);
  setEigenSparseVectorZero(m_NewtonInc);

  m_A.reserve(entriesPerCol);
  for (auto Constr : m_constraints) {
    m_A += Constr->getAMatrix(pointers);
  }

  m_A.makeCompressed();

  m_ATranspose = m_A.transpose();
}

auto ConstraintHandler::getAMatrix() -> Types::SparseMatrix<prec, indexType> & {
  return m_A;
}

auto ConstraintHandler::getATranspose()
    -> Types::SparseMatrix<prec, indexType> & {
  return m_ATranspose;
}

auto ConstraintHandler::getNewConstraint(ConstraintTypes constraintType)
    -> std::shared_ptr<BaseConstraint> {

  indexType pos = m_constraints.size();
  switch (constraintType) {
  case ConstraintTypes::GeneralLink:
    m_constraints.emplace_back();
    m_constraints[pos] = std::make_shared<GeneralLink>();
    m_constraints[pos]->setId(pos);
    break;
  default:
    break;
  }
  return m_constraints[pos];
}

auto ConstraintHandler::getConstraint(indexType constraintNumber)
    -> std::shared_ptr<BaseConstraint> {
  return m_constraints[constraintNumber];
}

auto ConstraintHandler::getNumberOfConstraints() -> indexType {
  return m_constraints.size();
}

void ConstraintHandler::modifyEquationSystem(
    PointerCollection &pointers, Eigen::SparseMatrix<prec, 0, indexType> &Kaa,
    Eigen::SparseMatrix<prec, 0, indexType> &Kab,
    Eigen::SparseMatrix<prec, 0, indexType> &Kba,
    Eigen::SparseMatrix<prec, 0, indexType> &Kbb, Types::VectorX<prec> &Fa,
    Types::VectorX<prec> &Fb, bool symmetric, bool upper) {
  if (m_constraints.size() != 0) {
    this->updateDB(pointers);
    if (symmetric) {
      if (upper) { // symmetric upper

        Kba += Kbb.selfadjointView<Eigen::Upper>() * m_A;
        Kaa += (m_ATranspose * Kba).triangularView<Eigen::Upper>();
        Kaa += (Kab * m_A).triangularView<Eigen::Upper>();
        Kab += m_ATranspose * Kbb.selfadjointView<Eigen::Upper>();

      }      // symmetric upper end
      else { // symmetric lower

        Kba += Kbb.selfadjointView<Eigen::Lower>() * m_A;
        Kaa += (m_ATranspose * Kba).triangularView<Eigen::Lower>();
        Kaa += (Kab * m_A).triangularView<Eigen::Lower>();
        Kab += m_ATranspose * Kbb.selfadjointView<Eigen::Lower>();

      }    // symmetric lower end
    }      // symmetric end
    else { // unsymmetric

      Kba += Kbb * m_A;
      Kaa += m_ATranspose * Kba;
      Kaa += Kab * m_A;
      Kab += m_ATranspose * Kbb;

    } // unsymmetric end
    Fa += m_ATranspose * Fb;
    Fa -= Kab * m_dB;
  }
}

auto ConstraintHandler::getSlaveNewtonSolution(PointerCollection &pointers)
    -> Types::SparseVector<prec, indexType> {
  setEigenSparseVectorZero(m_NewtonInc);

  for (auto i : m_constraints) {
    m_NewtonInc += i->getSlaveNewtonSolution(pointers);
  }

  return m_NewtonInc;
}

void ConstraintHandler::GeneralLinkGeo(
    PointerCollection &pointers, Geometry::GeometryTypes geoType,
    const std::vector<indexType> &masterNumbers,
    const std::vector<indexType> &slaveNumbers, indexType meshId,
    indexType order, indexType masterDof, indexType slaveDof, prec factor,
    prec difference, bool reorient) {
  if (masterNumbers.size() != slaveNumbers.size()) {
    pointers.getSPDLogger().error(
        "ConstraintHandler::GeneralLinkGeo: Number of masterelements does no "
        "math number of slaveelements");
    return;
  } else {
    indexType numGeo = masterNumbers.size();
    if (geoType == Geometry::GeometryTypes::Vertex) {
      for (indexType i = 0; i < numGeo; ++i) {
        this->LinkVertex(pointers, meshId, masterNumbers[i], slaveNumbers[i],
                         masterDof, slaveDof, factor, difference);
      }
    } else if (geoType == Geometry::GeometryTypes::Edges) {
      for (indexType i = 0; i < numGeo; ++i) {
        this->LinkEdges(pointers, meshId, masterNumbers[i], slaveNumbers[i],
                        masterDof, slaveDof, factor, difference, reorient);
      }
    } else if (geoType == Geometry::GeometryTypes::Faces) {
      for (indexType i = 0; i < numGeo; ++i) {
        this->LinkFaces(pointers, meshId, masterNumbers[i], slaveNumbers[i],
                        masterDof, slaveDof, factor, difference, reorient);
      }
    }
  }
}

void ConstraintHandler::toFile(std::ofstream &out) {
  writeEigenSparseMatrix(out, m_A);
  writeEigenSparseMatrix(out, m_ATranspose);
  writeEigenSparseVector(out, m_NewtonInc);
  writeEigenSparseVector(out, m_dB);

  indexType numConstr = m_constraints.size();
  writeScalar(out, numConstr);
  for (auto i : m_constraints) {
    writeScalar(out, i->getType());
    i->toFile(out);
  }
}

void ConstraintHandler::fromFile(std::ifstream &in) {
  readEigenSparseMatrix(in, m_A);
  readEigenSparseMatrix(in, m_ATranspose);
  readEigenSparseVector(in, m_NewtonInc);
  readEigenSparseVector(in, m_dB);

  indexType numConstr;
  readScalar(in, numConstr);
  m_constraints.reserve(numConstr);
  for (auto i = 0; i < numConstr; ++i) {
    ConstraintTypes coType;
    readScalar(in, coType);
    auto Constr = this->getNewConstraint(coType);
    Constr->fromFile(in);
  }
}

void ConstraintHandler::updateDB(PointerCollection &pointers) {
  setEigenSparseVectorZero(m_dB);

  for (auto i : m_constraints) {
    auto temp = i->getDB(pointers);
    m_dB += temp;
  }
}

void ConstraintHandler::LinkDofs(PointerCollection &pointers,
                                 indexType masterDof, indexType slaveDof,
                                 prec factor, prec difference) {
  auto &mDof = pointers.getEquationHandler()->getDegreeOfFreedom(masterDof);
  auto &sDof = pointers.getEquationHandler()->getDegreeOfFreedom(slaveDof);

  if (mDof.getStatus() == dofStatus::active &&
      sDof.getStatus() == dofStatus::active) {
    auto Constr = this->getNewConstraint(ConstraintTypes::GeneralLink);
    auto genLink = reinterpret_cast<std::shared_ptr<GeneralLink> &>(Constr);
    genLink->set(pointers, masterDof, slaveDof, factor, difference);
  }
}

void ConstraintHandler::LinkVertex(PointerCollection &pointers,
                                   indexType meshId, indexType masterVert,
                                   indexType slaveVert, indexType masterDof,
                                   indexType slaveDof, prec factor,
                                   prec difference) {

  auto &mVert = pointers.getGeometryData()->getVertexData(masterVert);
  auto &sVert = pointers.getGeometryData()->getVertexData(slaveVert);
  std::vector<GenericNodes *> mNodes, sNodes;
  mVert.getNodes(mNodes, meshId);
  sVert.getNodes(sNodes, meshId);

  this->LinkDofs(pointers, mNodes[0]->getDegreeOfFreedom(masterDof).getId(),
                 sNodes[0]->getDegreeOfFreedom(slaveDof).getId(), factor,
                 difference);
}

void ConstraintHandler::LinkEdges(PointerCollection &pointers, indexType meshId,
                                  indexType masterEdge, indexType slaveEdge,
                                  indexType masterDof, indexType slaveDof,
                                  prec factor, prec difference,
                                  bool reorient) {
  auto geoData = pointers.getGeometryData();
  auto &mEdge = geoData->getEdgeData(masterEdge);
  auto &sEdge = geoData->getEdgeData(slaveEdge);

  auto nvertsm = mEdge.getNumberOfVerts();

  if (mEdge.getType() != sEdge.getType()) {
    pointers.getSPDLogger().error(
        "Error when trying to link edge: {} to edge: {}. Edge types are not "
        "conforming! Aborting!!",
        slaveEdge, masterEdge);
    return;
  }

  IntegrationPoint ip;
  ip.xi = prec(0);
  auto mEvec = mEdge.getA1Vector(ip);
  auto sEvec = sEdge.getA1Vector(ip);
  prec dotVal = mEvec.dot(sEvec);
  if (dotVal > prec(1.01) || dotVal < prec(0.99)) {
    if (reorient) {
      sEdge.flip();
      geoData
          ->getEdgeRuntime(sEdge.getId())
          ->flip();
    } else {
      pointers.getSPDLogger().warn("Warning: When linking edge: {} to edge: {}."
                                   " Edge orientations are different! Check if "
                                   "this is intended.",
                                   slaveEdge, masterEdge);
    }
  }

  for (indexType i = 0; i < nvertsm; ++i) {
    this->LinkVertex(pointers, meshId, mEdge.getVertexNumber(i),
                     sEdge.getVertexNumber(i), masterDof, slaveDof, factor,
                     difference);
  }

  auto mNodes = mEdge.getNodesOfSet(meshId);
  auto sNodes = sEdge.getNodesOfSet(meshId);

  if (mNodes.size() == sNodes.size() && mNodes.size() != 0) {
    for (indexType i = 0; i < static_cast<indexType>(mNodes.size()); ++i) {
      auto &mDof = mNodes[i]->getDegreeOfFreedom(masterDof);
      auto &sDof = sNodes[i]->getDegreeOfFreedom(slaveDof);
      this->LinkDofs(pointers, mDof.getId(), sDof.getId(), factor, difference);
    }
  }
}

void ConstraintHandler::LinkFaces(PointerCollection &pointers, indexType meshId,
                                  indexType masterFace, indexType slaveFace,
                                  indexType masterDof, indexType slaveDof,
                                  prec factor, prec difference,
                                  bool reorient) {
  auto geoData = pointers.getGeometryData();
  auto mFace = geoData->getFaceData(masterFace);
  auto sFace = geoData->getFaceData(slaveFace);

  if (mFace->getType() != sFace->getType()) {
    pointers.getSPDLogger().error(
        "Error when trying to link face: {} to face: {}. Face types "
        "are not conforming! Aborting!!",
        slaveFace, masterFace);
    return;
  }

  if (reorient) {
    Types::Vector3<prec> snormal = sFace->getFaceNormal();
    Types::Vector3<prec> mnormal = mFace->getFaceNormal();
    if (abs(snormal.dot(mnormal) - prec(1)) >
        std::numeric_limits<prec>::epsilon() * prec(1000)) {
      sFace->flip();
      geoData->getFaceRuntime(sFace->getId())->flip();
    }
      
      

    IntegrationPoint faceIp;
    faceIp.xi = prec(0);
    faceIp.eta = prec(0);
    auto mG1 = mFace->getTangent_G1(faceIp).normalized();
    auto sG1 = sFace->getTangent_G1(faceIp).normalized();

    prec dd = abs(mG1.dot(sG1) - prec(1));
    indexType cc = 0;
    indexType numVerts = sFace->getNumberOfVerts();
    while (dd > std::numeric_limits<prec>::epsilon() * prec(1000) &&
           cc < numVerts) {
      sFace->rotate(1);
      geoData->getFaceRuntime(sFace->getId())->rotate(1);
      auto sG1 = sFace->getTangent_G1(faceIp).normalized();
      dd = abs(mG1.dot(sG1) - prec(1));
      ++cc;
    }
    if (dd > std::numeric_limits<prec>::epsilon() * prec(1000)) {
      pointers.getSPDLogger().warn(
          "When linking face: {} to face: {}. Face orientations are different! "
          "Make sure that everything is correct!!",
          slaveFace, masterFace);
    }
  }

  indexType nEdges = mFace->getNumberOfEdges();

  for (indexType i = 0; i < nEdges; ++i) {
    this->LinkEdges(pointers, meshId, mFace->getEdgeNumber(i), sFace->getEdgeNumber(i),
                    masterDof, slaveDof, factor, difference, reorient);
  }

  auto mNodes = mFace->getNodesOfSet(meshId);
  auto sNodes = sFace->getNodesOfSet(meshId);

  if (mNodes.size() == sNodes.size() && mNodes.size() != 0) {
    for (indexType i = 0; i < static_cast<indexType>(mNodes.size()); ++i) {
      auto &mDof = mNodes[i]->getDegreeOfFreedom(masterDof);
      auto &sDof = sNodes[i]->getDegreeOfFreedom(slaveDof);
      this->LinkDofs(pointers, mDof.getId(), sDof.getId(), factor, difference);
    }
  }
}

} // namespace HierAMuS