// Copyright 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "HomogenizationBeam.h"
#include "MatrixTypes.h"
#include "datatypes.h"
#include "solver/Constraints/GeneralLink.h"
#include "solver/Homogenization/HomogenizationBeam.h"

#include "equations/EquationHandler.h"
#include "geometry/GeometryData.h"
#include "geometry/GeometryTypes.h"

#include "control/BinaryWrite.h"

#include "spdlog/fmt/ostr.h"

namespace HierAMuS {
HomogenizationBeam::HomogenizationBeam()
    : xmin({0.0, 0.0, 0.0}), xmax({0.0, 0.0, 0.0}) {}

HomogenizationBeam::HomogenizationBeam(const HomogenizationBeam &other) : HomogenizationBase(other)
{
  this->homogenizationMatrix = other.homogenizationMatrix;
  this->meshIdDisp = other.meshIdDisp;
  this->dispOrder = other.dispOrder;
  this->bctype = other.bctype;
  this->leftFacesMaster = other.leftFacesMaster;
  this->rightFacesSlave = other.rightFacesSlave;
  this->xmin = other.xmin;
  this->xmax = other.xmax;
}

void HomogenizationBeam::init(PointerCollection &pointers,
                                          ParameterList &parameters) {

  meshIdDisp = parameters.getIndexVal("meshIdDisp");
  dispOrder = parameters.getIndexVal("dispOrder");
  bctype = parameters.getIndexVal("bctype");

  auto &Logger = pointers.getSPDLogger();
  Logger.info("Initializing HomogenizationBeam.....");
  Logger.info("    Mesh Id for displacements:                   {:>4}",meshIdDisp);
  Logger.info("    Order of approximation of displacements:     {:>4}",dispOrder);

  if (bctype == 0) {
    Logger.info("    Using displacement boundary conditions.");
  } else if (bctype == 1) {
    Logger.info("    Using periodic boundary conditions.");
  } else {
    Logger.warn("   Specified unsupported boundary condition type, using displacement boundary conditions.");
    bctype = 0;
  }

  auto getFaceNumbers = [](PointerCollection &pointers,
                           Types::Vector3<prec> &normal,
                           Types::Vector3<prec> &point) {
    auto faces =
        pointers.getGeometryData()->getFacesInPlane(pointers, normal, point);
    std::vector<indexType> faceNums;
    indexType nfaces = faces.size();
    faceNums.resize(nfaces);
    for (indexType i = 0; i < nfaces; ++i) {
      faceNums[i] = faces[i]->getId();
    }
    return faceNums;
  };
  this->xmax = pointers.getGeometryData()->getxMax();
  this->xmin = pointers.getGeometryData()->getxMin();

  Types::Vector3<prec> normal = {1.0, 0.0, 0.0};

  this->leftFacesMaster = getFaceNumbers(pointers, normal, this->xmin);
  this->rightFacesSlave = getFaceNumbers(pointers, normal, this->xmax);
  

  if (bctype == 0) {
    this->setDisplacementBoundaryConditions(pointers);
  } else if (bctype == 3) {
    this->setDisplacementBoundaryConditions2(pointers);
  } else if (bctype == 1) {
    pointers.getGeometryData()->sortReorientFacesPeriodicBC(
        pointers, this->leftFacesMaster, this->rightFacesSlave);
    auto &Logger = pointers.getSPDLogger();
    Logger.debug("Arranged master and slave faces for periodic boundary conditions");
    Logger.debug("Master faces: {}", fmt::join(leftFacesMaster, " "));
    Logger.debug("Master faces: {}", fmt::join(rightFacesSlave, " "));
    this->setPeriodicBoundaryConditions(pointers);
  }
}

void HomogenizationBeam::setDisplacementBoundaryConditions(
    PointerCollection &pointers) {

  Types::Vector3<indexType> dofs = {1, 1, 1};
  for (const auto &v :
       {std::cref(leftFacesMaster), std::cref(rightFacesSlave)}) {
    for (auto i : v.get()) {
      auto face = pointers.getGeometryData()->getFace(i);
      face->setBoundaryCondition(pointers, meshIdDisp, dispOrder,
                                 Geometry::ShapeFunctionTypes::H1, dofs, true);
    }
  }

}

void HomogenizationBeam::setDisplacementBoundaryConditions2(
    PointerCollection &pointers) {
  Types::Vector3<indexType> dofs = {1, 1, 1};
  for (const auto &v :
       {std::cref(leftFacesMaster), std::cref(rightFacesSlave)}) {
    for (auto i : v.get()) {
      auto face = pointers.getGeometryData()->getFace(i);
      face->setBoundaryCondition(pointers, meshIdDisp, dispOrder,
                                 Geometry::ShapeFunctionTypes::H1, dofs, true);
    }
  }
}

void HomogenizationBeam::setPeriodicBoundaryConditions(
    PointerCollection &pointers) {
  if (bctype == 1) {
    pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
        pointers, Geometry::GeometryTypes::Faces, this->leftFacesMaster,
        this->rightFacesSlave, this->meshIdDisp, this->dispOrder, 0, 0, 1, 0,
        true);
    pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
        pointers, Geometry::GeometryTypes::Faces, this->leftFacesMaster,
        this->rightFacesSlave, this->meshIdDisp, this->dispOrder, 1, 1, 1, 0,
        true);
    pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
        pointers, Geometry::GeometryTypes::Faces, this->leftFacesMaster,
        this->rightFacesSlave, this->meshIdDisp, this->dispOrder, 2, 2, 1, 0,
        true);

  }
}

void HomogenizationBeam::computeAMatrix(
    PointerCollection &pointers) {
  indexType totalEqs =
      pointers.getEquationHandler()->getNumberOfTotalEquations();
  indexType activeEqs =
      pointers.getEquationHandler()->getNumberOfActiveEquations();
  indexType inactiveIds = totalEqs - activeEqs;

  this->homogenizationMatrix.resize(inactiveIds, 6);
  this->homogenizationMatrix.setZero();

  if (bctype == 0) {
    for (const auto &ff :
         {std::cref(leftFacesMaster), std::cref(rightFacesSlave)}) {
      for (auto fNum : ff.get()) {
        auto face = pointers.getGeometryData()->getFace(fNum);
        indexType numV = face->getNumberOfVerts();
        for (indexType i = 0; i < numV; ++i) {
          auto vert = face->getVertex(pointers, i);
          auto Nodes = vert->getNodesOfSet(pointers, meshIdDisp);
          auto coor = vert->getCoordinates();
          for (auto node : Nodes) {
            indexType posA = node->getDegreeOfFreedom(0).getEqId();
            indexType posB = node->getDegreeOfFreedom(1).getEqId();
            indexType posC = node->getDegreeOfFreedom(2).getEqId();
            this->homogenizationMatrix(posA, 0) = coor(0); // x
            this->homogenizationMatrix(posB, 1) = coor(0); // x
            this->homogenizationMatrix(posC, 2) = coor(0); // x

            this->homogenizationMatrix(posB, 3) = -coor(0) * coor(2); // -x*z
            this->homogenizationMatrix(posC, 3) = coor(0) * coor(1);  // x*y
            this->homogenizationMatrix(posA, 4) = coor(0) * coor(2);  // x*z
            this->homogenizationMatrix(posA, 5) = -coor(0) * coor(1); // -x*y
          }
          
        }
      }
    }

  } else if (bctype == 3) {
    for (const auto &ff :
         {std::cref(leftFacesMaster), std::cref(rightFacesSlave)}) {
      for (auto fNum : ff.get()) {
        auto face = pointers.getGeometryData()->getFace(fNum);
        indexType numV = face->getNumberOfVerts();
        for (indexType i = 0; i < numV; ++i) {
          auto vert = face->getVertex(pointers, i);
          auto Nodes = vert->getNodesOfSet(pointers, meshIdDisp);
          auto coor = vert->getCoordinates();
          for (auto node : Nodes) {
            indexType posA = node->getDegreeOfFreedom(0).getEqId();
            indexType posB = node->getDegreeOfFreedom(1).getEqId();
            indexType posC = node->getDegreeOfFreedom(2).getEqId();
            this->homogenizationMatrix(posA, 0) = coor(0); // x
            this->homogenizationMatrix(posB, 1) = coor(0); // x
            this->homogenizationMatrix(posC, 2) = coor(0); // x

            this->homogenizationMatrix(posB, 3) = -coor(0) * coor(2); // -x*z
            this->homogenizationMatrix(posC, 3) = coor(0) * coor(1);  // x*y
            this->homogenizationMatrix(posA, 4) = coor(0) * coor(2);  // x*z
            this->homogenizationMatrix(posA, 5) = -coor(0) * coor(1); // -x*y
          }
        }
      }
    }

  } else if (bctype == 1) {
    for (indexType i = 0; i < leftFacesMaster.size(); ++i) {
      auto mFace = pointers.getGeometryData()->getFace(leftFacesMaster[i]);
      auto sFace = pointers.getGeometryData()->getFace(rightFacesSlave[i]);
      indexType numv = mFace->getNumberOfVerts();
      for (indexType lv = 0; lv < numv; ++lv) {
        auto mVert = mFace->getVertex(pointers, lv);
        auto sVert = sFace->getVertex(pointers, lv);
        Types::Vector3<prec> dx =
            sVert->getCoordinates() - mVert->getCoordinates();

        Types::Vector3<prec> coor = sVert->getCoordinates();
        auto sNodes = sVert->getNodesOfSet(pointers, meshIdDisp);
        indexType posA = sNodes[0]->getDegreeOfFreedom(0).getEqId();
        indexType posB = sNodes[0]->getDegreeOfFreedom(1).getEqId();
        indexType posC = sNodes[0]->getDegreeOfFreedom(2).getEqId();
        this->homogenizationMatrix(posA, 0) = dx(0); // x
        this->homogenizationMatrix(posB, 1) = dx(0); // x
        this->homogenizationMatrix(posC, 2) = dx(0); // x

        this->homogenizationMatrix(posB, 3) = -dx(0) * coor(2); // -x*z
        this->homogenizationMatrix(posC, 3) = dx(0) * coor(1);  // x*y
        this->homogenizationMatrix(posA, 4) = dx(0) * coor(2);  // x*z
        this->homogenizationMatrix(posA, 5) = -dx(0) * coor(1); // -x*y
      }
    }
  }

  pointers.getSPDLogger().trace("Computed A-Matrix:\n{}",this->homogenizationMatrix);
}

auto HomogenizationBeam::getNumberOfStrains() -> indexType {
  return indexType(6);
}

auto HomogenizationBeam::getAMatrix() -> Types::MatrixXX<prec> & {
  return this->homogenizationMatrix;
}

auto HomogenizationBeam::getDv() -> prec {
  return (xmax(0) - xmin(0));
}

auto HomogenizationBeam::getDisplacementIncrement(
    Types::VectorX<prec> &strainIncrement) -> Types::VectorX<prec> {
  Types::VectorX<prec> disp;
  if (bctype == 0) {
    disp = this->homogenizationMatrix * strainIncrement;
  } else if (bctype == 1) {
    Types::VectorX<prec> sttemp = strainIncrement;
    for (indexType i = 0; i < 6; ++i) {
      sttemp(i) = prec(0);
    }
    disp = homogenizationMatrix * sttemp;
  }

  return disp;
}

void HomogenizationBeam::setPeriodicDisplacements(
    PointerCollection &pointers, Types::VectorX<prec> &strains,
    Types::VectorX<prec> &strainIncrement) {

  if (bctype == 1) {
    auto setDb = [](PointerCollection &pointers, DegreeOfFreedom &Dof, prec B) {
      if (Dof.getStatus() == dofStatus::constraint) {
        auto bConst =
            pointers.getSolutionState()->getConstraint(Dof.getConstraintId());
        auto genLink = reinterpret_cast<std::shared_ptr<GeneralLink> &>(bConst);
        genLink->setB(B);
      }
    };
    Types::Vector3<prec> u;
    Types::Matrix3X<prec> A;
    A.resize(3, strains.rows());
    A.setZero();
    for (indexType i = 0; i < leftFacesMaster.size(); ++i) {
      auto mFace = pointers.getGeometryData()->getFace(leftFacesMaster[i]);
      auto sFace = pointers.getGeometryData()->getFace(rightFacesSlave[i]);

      indexType numV = mFace->getNumberOfVerts();
      for (indexType nn = 0; nn < numV; ++nn) {
        auto mVert = mFace->getVertex(pointers, nn);
        auto sVert = sFace->getVertex(pointers, nn);
        Types::Vector3<prec> sCoor = sVert->getCoordinates();
        Types::Vector3<prec> dx = sCoor - mVert->getCoordinates();
        auto sNodes = sVert->getNodesOfSet(pointers, meshIdDisp);
        A(0, 0) = dx(0);
        A(1, 1) = dx(0);
        A(2, 2) = dx(0);
        A(1, 3) = -dx(0) * sCoor(2); // -x*z
        A(2, 3) = dx(0) * sCoor(1);  // x*y
        A(0, 4) = dx(0) * sCoor(2);  // x*z
        A(0, 5) = -dx(0) * sCoor(1); // -x*y
        u = A * strains;
        auto &DofA = sNodes[0]->getDegreeOfFreedom(0);
        auto &DofB = sNodes[0]->getDegreeOfFreedom(1);
        auto &DofC = sNodes[0]->getDegreeOfFreedom(2);
        setDb(pointers, DofA, u(0));
        setDb(pointers, DofB, u(1));
        setDb(pointers, DofC, u(2));
      }
    }
  }
}

void HomogenizationBeam::toFile(PointerCollection &pointers,
                                            std::ofstream &out) {
  HomogenizationBase::toFile(pointers, out);
  writeScalar(out, meshIdDisp);
  writeScalar(out, dispOrder);
  writeScalar(out, bctype);
  writeStdVector(out, leftFacesMaster);
  writeStdVector(out, rightFacesSlave);

  writeEigenMatrix(out, homogenizationMatrix);
  writeEigenMatrix(out, xmin);
  writeEigenMatrix(out, xmax);
  
}

void HomogenizationBeam::fromFile(PointerCollection &pointers,
                                              std::ifstream &in) {
  HomogenizationBase::fromFile(pointers, in);
  readScalar(in, meshIdDisp);
  readScalar(in, dispOrder);
  readScalar(in, bctype);
  readStdVector(in, leftFacesMaster);
  readStdVector(in, rightFacesSlave);

  readEigenMatrix(in, homogenizationMatrix);
  readEigenMatrix(in, xmin);
  readEigenMatrix(in, xmax);
}

} // namespace HierAMuS