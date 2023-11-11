
#include "HomogenizationShell.h"
#include "MatrixTypes.h"
#include "datatypes.h"
#include "solver/Homogenization/HomogenizationShell.h"

//Equations
#include "EquationHandler.h"

#include "geometry/GeometryData.h"
#include "geometry/GeometryTypes.h"
#include "geometry/Faces/FacesData.h"

#include "control/BinaryWrite.h"

#include "spdlog/fmt/ostr.h"

namespace HierAMuS {
HomogenizationShell::HomogenizationShell()
    : xmin({0.0, 0.0, 0.0}), xmax({0.0, 0.0, 0.0}) {}

HomogenizationShell::HomogenizationShell(const HomogenizationShell &other) : HomogenizationBase(other)
{
  this->homogenizationMatrix = other.homogenizationMatrix;

  this->meshIdDisp = other.meshIdDisp;
  this->dispOrder = other.dispOrder;
  this->bctype = other.bctype;

  this->leftFaces = other.leftFaces;
  this->rightFaces = other.rightFaces;
  this->frontFaces = other.frontFaces;
  this->backFaces = other.backFaces;

  this->xmin = other.xmin;
  this->xmax = other.xmax;

}

void HomogenizationShell::init(PointerCollection &pointers,
                                 ParameterList &parameters) {


  meshIdDisp = parameters.getIndexVal("meshIdDisp");
  dispOrder = parameters.getIndexVal("dispOrder");
  bctype = parameters.getIndexVal("bctype");

  auto &Logger = pointers.getSPDLogger();
  Logger.info("Initializing HomogenizationShell.....");
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

  
  auto getFaceNumbers = [](PointerCollection &pointers, Types::Vector3<prec> &normal, Types::Vector3<prec> & point)
  {
    auto faces =
        pointers.getGeometryData()->getFacesInPlane(normal, point);
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
  
  this->leftFaces = getFaceNumbers(pointers,normal,this->xmin);
  this->rightFaces = getFaceNumbers(pointers,normal,this->xmax);
  
  Types::Vector3<prec> normal1 = {0.0, 1.0, 0.0};
  
  this->backFaces = getFaceNumbers(pointers,normal1,this->xmin);
  this->frontFaces = getFaceNumbers(pointers,normal1,this->xmax);

  this->setDisplacementBoundaryConditions(pointers);
  
}



void HomogenizationShell::setDisplacementBoundaryConditions(
    PointerCollection &pointers) {

  Types::Vector3<indexType> dofs = {1, 1, 1};
  for(auto i : leftFaces){
    auto face = pointers.getGeometryData()->getFaceData(i);
    face->setBoundaryCondition(meshIdDisp, dispOrder, Geometry::ShapeFunctionTypes::H1, dofs, true);
  }
  for(auto i : rightFaces){
    auto face = pointers.getGeometryData()->getFaceData(i);
    face->setBoundaryCondition(meshIdDisp, dispOrder, Geometry::ShapeFunctionTypes::H1, dofs, true);
  }
  for(auto i : backFaces){
    auto face = pointers.getGeometryData()->getFaceData(i);
    face->setBoundaryCondition(meshIdDisp, dispOrder, Geometry::ShapeFunctionTypes::H1, dofs, true);
  }
  for(auto i : frontFaces){
    auto face = pointers.getGeometryData()->getFaceData(i);
    face->setBoundaryCondition(meshIdDisp, dispOrder, Geometry::ShapeFunctionTypes::H1, dofs, true);
  }

}

void HomogenizationShell::computeAMatrix(PointerCollection &pointers) {
  indexType totalEqs = pointers.getEquationHandler()->getNumberOfTotalEquations();
  indexType activeEqs = pointers.getEquationHandler()->getNumberOfActiveEquations();
  indexType inactiveIds = totalEqs - activeEqs;

  this->homogenizationMatrix.resize(inactiveIds, 8);//check
  this->homogenizationMatrix.setZero();
  std::vector<std::vector<indexType>> ffaces = { leftFaces, rightFaces, backFaces, frontFaces };
  if (bctype == 0) {
      for (auto j = 0; j < ffaces.size(); j++) {
          for (auto fNum : ffaces[j]) {
              auto face = pointers.getGeometryData()->getFaceData(fNum);
              indexType numV = face->getNumberOfVerts();
              for (indexType i = 0; i < numV; ++i) {
                  auto vert = face->getVertex(i);
                  auto Nodes = vert->getNodesOfSet(meshIdDisp);
                  auto coor = vert->getCoordinates();
                  for (auto node : Nodes) {

                      indexType posA = node->getDegreeOfFreedom(0).getEqId();
                      indexType posB = node->getDegreeOfFreedom(1).getEqId();

                      this->homogenizationMatrix(posA, 0) = coor(0);                           //xx
                      this->homogenizationMatrix(posA, 2) = coor(1) * prec(0.5);               //0.5y
                      this->homogenizationMatrix(posA, 3) = coor(0) * coor(2);                 //xz
                      this->homogenizationMatrix(posA, 5) = coor(1) * coor(2) * prec(0.5);     //0.5yz
                      this->homogenizationMatrix(posA, 6) = coor(2);                           //z

                      this->homogenizationMatrix(posB, 1) = coor(1);                           //yy
                      this->homogenizationMatrix(posB, 2) = coor(0) * prec(0.5);               //0.5x
                      this->homogenizationMatrix(posB, 4) = coor(1) * coor(2);                 //yz
                      this->homogenizationMatrix(posB, 5) = coor(0) * coor(2) * prec(0.5);     //0.5xz
                      this->homogenizationMatrix(posB, 7) = coor(2);                           //z

                  }
              }
          }
      }
  }
    
  pointers.getSPDLogger().trace("Computed A-Matrix:\n{}",this->homogenizationMatrix);

}

auto HomogenizationShell::getNumberOfStrains() -> indexType {
  return indexType(8);
}

auto HomogenizationShell::getAMatrix() -> Types::MatrixXX<prec> & {
  return this->homogenizationMatrix;
}

auto HomogenizationShell::getDv() -> prec { return (xmax(0) - xmin(0))*(xmax(1) - xmin(1));}



auto HomogenizationShell::getDisplacementIncrement(
    Types::VectorX<prec> &strainIncrement) -> Types::VectorX<prec> {
  Types::VectorX<prec> disp =
      this->homogenizationMatrix * strainIncrement;
  return disp;
}

void HomogenizationShell::toFile(PointerCollection &pointers,
                                   std::ofstream &out) {
  HomogenizationBase::toFile(pointers, out);
  writeScalar(out, meshIdDisp);
  writeScalar(out, dispOrder);
  writeScalar(out, bctype);
  writeStdVector(out, leftFaces);
  writeStdVector(out, rightFaces);
   writeStdVector(out, backFaces);
  writeStdVector(out, frontFaces);

  writeEigenMatrix(out, homogenizationMatrix);
  writeEigenMatrix(out,xmin);
  writeEigenMatrix(out,xmax);

}

void HomogenizationShell::fromFile(PointerCollection &pointers,
                                     std::ifstream &in) {
  HomogenizationBase::fromFile(pointers, in);
  readScalar(in, meshIdDisp);
  readScalar(in, dispOrder);
  readScalar(in, bctype);
  readStdVector(in, leftFaces);
  readStdVector(in, rightFaces);
  readStdVector(in, backFaces);
  readStdVector(in, frontFaces);

  readEigenMatrix(in, homogenizationMatrix);
  readEigenMatrix(in,xmin);
  readEigenMatrix(in,xmax);
}

} // namespace FEMProject