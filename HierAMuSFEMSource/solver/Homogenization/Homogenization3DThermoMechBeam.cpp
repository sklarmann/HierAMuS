// Copyright 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "Homogenization3DThermoMechBeam.h"
#include "MatrixTypes.h"
#include "datatypes.h"
#include "solver/Homogenization/Homogenization3DThermoMechBeam.h"
#include "solver/Constraints/GeneralLink.h"

#include "equations/EquationHandler.h"
#include "geometry/GeometryData.h"
#include "geometry/GeometryTypes.h"

#include "control/BinaryWrite.h"

#include "spdlog/fmt/ostr.h"
#include <spdlog/fmt/bundled/format.h>

namespace HierAMuS {
Homogenization3DThermoMechBeam::Homogenization3DThermoMechBeam()
    : xmin({0.0, 0.0, 0.0}), xmax({0.0, 0.0, 0.0}) {}

Homogenization3DThermoMechBeam::Homogenization3DThermoMechBeam(
    const Homogenization3DThermoMechBeam &other) : HomogenizationBase(other)
{
  this->homogenizationMatrix = other.homogenizationMatrix;
  this->meshIdDisp = other.meshIdDisp;
  this->dispOrder = other.dispOrder;
  this->meshIdTemp = other.meshIdTemp;
  this->tempOrder = other.tempOrder;
  this->bctype = other.bctype;

  this->leftFacesMaster = other.leftFacesMaster;
  this->rightFacesSlave = other.rightFacesSlave;
  this->bottomFacesMaster = other.bottomFacesMaster;
  this->topFacesSlave = other.topFacesSlave;
  this->backFacesMaster = other.backFacesMaster;
  this->frontFacesSlave = other.frontFacesSlave;

  this->xmin = other.xmin;
  this->xmax = other.xmax;
}

void Homogenization3DThermoMechBeam::init(PointerCollection &pointers,
                                 ParameterList &parameters) {


  meshIdDisp = parameters.getIndexVal("meshIdDisp");
  dispOrder = parameters.getIndexVal("dispOrder");
  meshIdTemp = parameters.getIndexVal("meshIdTemp");
  tempOrder = parameters.getIndexVal("tempOrder");
  bctype = parameters.getIndexVal("bctype");

  auto Logger = pointers.getSPDLogger();

  

  if (bctype == 0) {
    Logger.info("\n{:-<100}\n"
              "Initializing Homogenization3DThermoMechBeam.....\n"
              "   Mesh Id for displacements:                 {:>12}\n"
              "   Order of approximation of displacements:   {:>12}\n"
              "   Using displacement boundary conditions."
              "\n{:-<100}\n",
              "",
              meshIdDisp,
              dispOrder,
              "");

  } else if (bctype == 1) {
    Logger.info("\n{:-<100}\n"
              "Initializing Homogenization3DThermoMechBeam.....\n"
              "   Mesh Id for displacements:                 {:>12}\n"
              "   Order of approximation of displacements:   {:>12}\n"
              "   Using periodic boundary conditions."
              "\n{:-<100}\n",
              "",
              meshIdDisp,
              dispOrder,
              "");
  } else {
    
    Logger.warn("\n{:-<100}\n"
              "Initializing Homogenization3DThermoMechBeam.....\n"
              "   Mesh Id for displacements:                 {:>12}\n"
              "   Order of approximation of displacements:   {:>12}\n"
              "   Specified unsupported boundary condition type, using displacement boundary conditions."
              "\n{:-<100}\n",
              "",
              meshIdDisp,
              dispOrder,
              "");
    bctype = 0;
  }

  this->computeGeometryParameters(pointers);

  if (bctype == 0) {
    this->setDisplacementBoundaryConditions(pointers);
  } else if (bctype == 3) 
  {
    this->setDisplacementBoundaryConditions2(pointers);
  }
  else if (bctype == 1) 
  {
    pointers.getGeometryData()->sortReorientFacesPeriodicBC(
        pointers, this->leftFacesMaster, this->rightFacesSlave);
    Logger.debug("Arranged master and slave faces for periodic boundary conditions\n"
                               "Master faces: {}\n"
                               "Slave faces: {}\n",fmt::join(leftFacesMaster, " "),fmt::join(rightFacesSlave, " "));
    this->setPeriodicBoundaryConditions(pointers);
  }
}

void Homogenization3DThermoMechBeam::computeGeometryParameters(PointerCollection &pointers){
  
  // Setting up the faces
  auto getFaceNumbers = [](PointerCollection &pointers, Types::Vector3<prec> &normal, Types::Vector3<prec> & point)
  {
    auto faces =
        pointers.getGeometryData()->getFacesInPlane(pointers, normal, point);
    std::vector<indexType> faceNums;
    indexType nfaces = faces.size();
    faceNums.resize(nfaces);
    for (indexType i=0;i<nfaces;++i)
    {
      faceNums[i] = faces[i]->getId();
    }
    return faceNums;
  };
  this->xmax = pointers.getGeometryData()->getxMax();
  this->xmin = pointers.getGeometryData()->getxMin();

  
  Types::Vector3<prec> normal = {1.0, 0.0, 0.0};

  this->leftFacesMaster = getFaceNumbers(pointers, normal, this->xmin);
  this->rightFacesSlave = getFaceNumbers(pointers, normal, this->xmax);

  normal = {0.0, 1.0, 0.0};
  this->backFacesMaster = getFaceNumbers(pointers, normal, this->xmin);
  this->frontFacesSlave = getFaceNumbers(pointers, normal, this->xmax);

  normal = {0.0, 0.0, 1.0};
  this->bottomFacesMaster = getFaceNumbers(pointers, normal, this->xmin);
  this->topFacesSlave = getFaceNumbers(pointers, normal, this->xmax);


  this->v1 = pointers.getGeometryData()->getVertexClosestTo({this->xmin[0], this->xmin[1], this->xmin[2]}).getId();
  this->v2 = pointers.getGeometryData()->getVertexClosestTo({this->xmax[0], this->xmin[1], this->xmin[2]}).getId();
  this->v3 = pointers.getGeometryData()->getVertexClosestTo({this->xmax[0], this->xmax[1], this->xmin[2]}).getId();
  this->v4 = pointers.getGeometryData()->getVertexClosestTo({this->xmin[0], this->xmax[1], this->xmin[2]}).getId();
  this->v5 = pointers.getGeometryData()->getVertexClosestTo({this->xmin[0], this->xmin[1], this->xmax[2]}).getId();
  this->v6 = pointers.getGeometryData()->getVertexClosestTo({this->xmax[0], this->xmin[1], this->xmax[2]}).getId();
  this->v7 = pointers.getGeometryData()->getVertexClosestTo({this->xmax[0], this->xmax[1], this->xmax[2]}).getId();
  this->v8 = pointers.getGeometryData()->getVertexClosestTo({this->xmin[0], this->xmax[1], this->xmax[2]}).getId();



}

void Homogenization3DThermoMechBeam::setDisplacementBoundaryConditions(
    PointerCollection &pointers) {

  Types::Vector3<indexType> dofs = {1, 1, 1};
  for (const auto &v : {std::cref(leftFacesMaster), std::cref(rightFacesSlave)}) {
    for (auto i : v.get()) {
      auto face = pointers.getGeometryData()->getFace(i);
      face->setBoundaryCondition(pointers, meshIdDisp, dispOrder,
                                 Geometry::ShapeFunctionTypes::H1, dofs, true);
      face->setBoundaryCondition(pointers, meshIdTemp, tempOrder,
                                 Geometry::ShapeFunctionTypes::H1, dofs, true);
    }
  }

  for (const auto &v :
       {std::cref(bottomFacesMaster), std::cref(topFacesSlave), std::cref(backFacesMaster), std::cref(frontFacesSlave)}) {
    for (auto i : v.get()) {
      auto face = pointers.getGeometryData()->getFace(i);
      face->setBoundaryCondition(pointers, meshIdTemp, tempOrder,
                                 Geometry::ShapeFunctionTypes::H1, dofs, true);
    }
  }

}

void Homogenization3DThermoMechBeam::setDisplacementBoundaryConditions2(
    PointerCollection &pointers)
{
  Types::Vector3<indexType> dofs = {1, 1, 1};
  for (const auto &v :
       {std::cref(leftFacesMaster), std::cref(rightFacesSlave)}) {
    for (auto i : v.get()) {
      auto face = pointers.getGeometryData()->getFace(i);
      face->setBoundaryCondition(pointers, meshIdDisp, dispOrder,
                                 Geometry::ShapeFunctionTypes::H1, dofs, true);
      face->setBoundaryCondition(pointers, meshIdTemp, tempOrder,
                                 Geometry::ShapeFunctionTypes::H1, dofs, true);
    }
  }

}

void Homogenization3DThermoMechBeam::setPeriodicBoundaryConditions(
    PointerCollection &pointers) {
  if(bctype == 1){
  pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
    pointers, Geometry::GeometryTypes::Faces, this->leftFacesMaster,
    this->rightFacesSlave, this->meshIdDisp, this->dispOrder, 0, 0, 1, 0, true);
  pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
    pointers, Geometry::GeometryTypes::Faces, this->leftFacesMaster,
    this->rightFacesSlave, this->meshIdDisp, this->dispOrder, 1, 1, 1, 0, true);
  pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
    pointers, Geometry::GeometryTypes::Faces, this->leftFacesMaster,
    this->rightFacesSlave, this->meshIdDisp, this->dispOrder, 2, 2, 1, 0, true);

  //Types::Vector3<indexType> dofs = {1, 1, 1};
  //for (const auto &v :
  //     {std::cref(bottomFacesMaster), std::cref(topFacesSlave),
  //      std::cref(backFacesMaster), std::cref(frontFacesSlave),
  //      std::cref(leftFacesMaster), std::cref(rightFacesSlave)}) {
  //  for (auto i : v.get()) {
  //    auto face = pointers.getGeometryData()->getFace(i);
  //    face->setBoundaryCondition(pointers, meshIdTemp, tempOrder,
  //                               Geometry::ShapeFunctionTypes::H1, dofs, true);
  //  }
  //}
  Types::Vector3<indexType> dofs = {1, 1, 1};
  for (const auto &v :
       {std::cref(leftFacesMaster), std::cref(rightFacesSlave),
        std::cref(topFacesSlave), std::cref(bottomFacesMaster),
        std::cref(backFacesMaster), std::cref(frontFacesSlave)}) {
    for (auto i : v.get()) {
      auto face = pointers.getGeometryData()->getFace(i);
      face->setBoundaryCondition(pointers, meshIdTemp, tempOrder,
                                 Geometry::ShapeFunctionTypes::H1, dofs, true);
    }
  }
  }
}

void Homogenization3DThermoMechBeam::computeAMatrix(PointerCollection &pointers) {
  indexType totalEqs = pointers.getEquationHandler()->getNumberOfTotalEquations();
  indexType activeEqs = pointers.getEquationHandler()->getNumberOfActiveEquations();
  indexType inactiveIds = totalEqs - activeEqs;

  this->homogenizationMatrix.resize(inactiveIds, 12);
  this->homogenizationMatrix.setZero();

  if(bctype==0){
  for (const auto &ff : {std::cref(leftFacesMaster), std::cref(rightFacesSlave)})
  {
    for (auto fNum : ff.get())
    {
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
        auto NodesT = vert->getNodesOfSet(pointers, meshIdTemp);
        for (auto node : NodesT) {
          indexType posA = node->getDegreeOfFreedom(0).getEqId();
          this->homogenizationMatrix(posA, 6) = prec(1.0);           // x
          this->homogenizationMatrix(posA, 7) = coor(2);             // z
          this->homogenizationMatrix(posA, 8) = -coor(1);            // -y
          this->homogenizationMatrix(posA, 9) = coor(0);             // x
          this->homogenizationMatrix(posA, 10) = coor(0) * coor(2);  // x*z
          this->homogenizationMatrix(posA, 11) = -coor(0) * coor(1); // -x*y
        }
      }
    }
  }
  for (const auto &ff :
       {std::cref(bottomFacesMaster), std::cref(topFacesSlave), std::cref(frontFacesSlave), std::cref(backFacesMaster)}) {
    for (auto fNum : ff.get()) {
      auto face = pointers.getGeometryData()->getFace(fNum);
      indexType numV = face->getNumberOfVerts();
      for (indexType i = 0; i < numV; ++i) {
        auto vert = face->getVertex(pointers, i);
        auto Nodes = vert->getNodesOfSet(pointers, meshIdTemp);
        auto coor = vert->getCoordinates();
        for (auto node : Nodes) {
          indexType posA = node->getDegreeOfFreedom(0).getEqId();
          this->homogenizationMatrix(posA, 6) = prec(1.0);           // 1
          this->homogenizationMatrix(posA, 7) = coor(2);             // z
          this->homogenizationMatrix(posA, 8) = -coor(1);            // -y
          this->homogenizationMatrix(posA, 9) = coor(0);             // x
          this->homogenizationMatrix(posA, 10) = coor(0) * coor(2);  // x*z
          this->homogenizationMatrix(posA, 11) = -coor(0) * coor(1); // -x*y
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
        auto NodesT = vert->getNodesOfSet(pointers, meshIdTemp);
        for (auto node : NodesT) {
          indexType posA = node->getDegreeOfFreedom(0).getEqId();
          this->homogenizationMatrix(posA, 6) = prec(1.0); // x
          this->homogenizationMatrix(posA, 7) = coor(2);             // x*z
          this->homogenizationMatrix(posA, 8) = -coor(1);            // -x*y
          this->homogenizationMatrix(posA, 9) = coor(0);             // x
          this->homogenizationMatrix(posA, 10) = coor(0) * coor(2);  // x*z
          this->homogenizationMatrix(posA, 11) = -coor(0) * coor(1); // -x*y
        }
      }
    }
  }
    
  }
  else if (bctype == 1) {
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
        this->homogenizationMatrix(posA, 0) = dx(0);   // x
        this->homogenizationMatrix(posB, 1) = dx(0);   // x
        this->homogenizationMatrix(posC, 2) = dx(0);   // x

        
        this->homogenizationMatrix(posB, 3) = -dx(0) * coor(2); // -x*z
        this->homogenizationMatrix(posC, 3) = dx(0) * coor(1);  // x*y
        this->homogenizationMatrix(posA, 4) = dx(0) * coor(2);  // x*z
        this->homogenizationMatrix(posA, 5) = -dx(0) * coor(1); // -x*y
      }
    }
    for (const auto &ff :
         {std::cref(leftFacesMaster), std::cref(rightFacesSlave)}) {
      for (auto fNum : ff.get()) {
        auto face = pointers.getGeometryData()->getFace(fNum);
        indexType numV = face->getNumberOfVerts();
        for (indexType i = 0; i < numV; ++i) {
        auto vert = face->getVertex(pointers, i);
        auto Nodes = vert->getNodesOfSet(pointers, meshIdTemp);
        auto coor = vert->getCoordinates();
        for (auto node : Nodes) {
          indexType posA = node->getDegreeOfFreedom(0).getEqId();
          this->homogenizationMatrix(posA, 6) = prec(1.0);           // x
          this->homogenizationMatrix(posA, 7) = coor(2);   // x*z
          this->homogenizationMatrix(posA, 8) = -coor(1);  // -x*y
          this->homogenizationMatrix(posA, 9) = coor(0);             // x
          this->homogenizationMatrix(posA, 10) = coor(0) * coor(2);  // x*z
          this->homogenizationMatrix(posA, 11) = -coor(0) * coor(1); // -x*y
        }
        }
      }
    }
    for (const auto &ff :
         {
          std::cref(topFacesSlave), std::cref(bottomFacesMaster),
          std::cref(frontFacesSlave), std::cref(backFacesMaster)}) {
      for (auto fNum : ff.get()) {
        auto face = pointers.getGeometryData()->getFace(fNum);
        indexType numV = face->getNumberOfVerts();
        for (indexType i = 0; i < numV; ++i) {
        auto vert = face->getVertex(pointers, i);
        auto Nodes = vert->getNodesOfSet(pointers, meshIdTemp);
        auto coor = vert->getCoordinates();
        for (auto node : Nodes) {
          indexType posA = node->getDegreeOfFreedom(0).getEqId();
          this->homogenizationMatrix(posA, 6) = prec(1.0);           // x
          this->homogenizationMatrix(posA, 7) = coor(2);   // x*z
          this->homogenizationMatrix(posA, 8) = -coor(1);  // -x*y
          this->homogenizationMatrix(posA, 9) = coor(0);             // x
          this->homogenizationMatrix(posA, 10) = coor(0) * coor(2);  // x*z
          this->homogenizationMatrix(posA, 11) = -coor(0) * coor(1); // -x*y
        }
        }
      }
    }
  }



  pointers.getSPDLogger().trace("Computed A-Matrix:\n{}",this->homogenizationMatrix);
}

auto Homogenization3DThermoMechBeam::getNumberOfStrains() -> indexType {
  return indexType(12);
}

auto Homogenization3DThermoMechBeam::getAMatrix() -> Types::MatrixXX<prec> & {
  return this->homogenizationMatrix;
}

auto Homogenization3DThermoMechBeam::getDv() -> prec { return (xmax(0) - xmin(0)); }



auto Homogenization3DThermoMechBeam::getDisplacementIncrement(
    Types::VectorX<prec> &strainIncrement) -> Types::VectorX<prec> {
  Types::VectorX<prec> disp;
  if (bctype == 0)
  {
    disp = this->homogenizationMatrix * strainIncrement;
  } else if (bctype == 1)
  {
    Types::VectorX<prec> sttemp = strainIncrement;
    for (indexType i=0;i<6;++i)
    {
      sttemp(i) = prec(0);
    }
    disp = homogenizationMatrix * sttemp;
  }
      
  return disp;
}

void Homogenization3DThermoMechBeam::setPeriodicDisplacements(
  PointerCollection& pointers, Types::VectorX<prec> &strains, Types::VectorX<prec> &strainIncrement)
{


  if (bctype==1)
  {
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
    for (indexType i=0;i<leftFacesMaster.size();++i) {
      auto mFace = pointers.getGeometryData()->getFace(leftFacesMaster[i]);
      auto sFace = pointers.getGeometryData()->getFace(rightFacesSlave[i]);

      indexType numV = mFace->getNumberOfVerts();
      for (indexType nn =0;nn<numV;++nn)
      {
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

void Homogenization3DThermoMechBeam::toFile(PointerCollection &pointers,
                                   std::ofstream &out) {
  HomogenizationBase::toFile(pointers, out);
  writeScalar(out, meshIdDisp);
  writeScalar(out, dispOrder);
  writeScalar(out, meshIdTemp);
  writeScalar(out, tempOrder);
  writeScalar(out, bctype);
  writeStdVector(out, leftFacesMaster);
  writeStdVector(out, rightFacesSlave);

  writeEigenMatrix(out, homogenizationMatrix);
  writeEigenMatrix(out,xmin);
  writeEigenMatrix(out,xmax);

}

void Homogenization3DThermoMechBeam::fromFile(PointerCollection &pointers,
                                     std::ifstream &in) {
  HomogenizationBase::fromFile(pointers, in);
  readScalar(in, meshIdDisp);
  readScalar(in, dispOrder);
  readScalar(in, meshIdTemp);
  readScalar(in, tempOrder);
  readScalar(in, bctype);
  readStdVector(in, leftFacesMaster);
  readStdVector(in, rightFacesSlave);

  readEigenMatrix(in, homogenizationMatrix);
  readEigenMatrix(in,xmin);
  readEigenMatrix(in,xmax);
}

} // namespace HierAMuS