
#include "Homogenization3DSolid.h"
#include "MatrixTypes.h"
#include "datatypes.h"
#include "solver/Homogenization/Homogenization3DSolid.h"

#include "equations/EquationHandler.h"
#include "geometry/GeometryData.h"
#include "geometry/GeometryTypes.h"

#include "control/BinaryWrite.h"

#include "spdlog/fmt/ostr.h"

namespace HierAMuS {
Homogenization3DSolid::Homogenization3DSolid()
    : xmin({0.0, 0.0, 0.0}), xmax({0.0, 0.0, 0.0}) {}

Homogenization3DSolid::Homogenization3DSolid(const Homogenization3DSolid &other)
    : HomogenizationBase(other) {
  this->homogenizationMatrix = other.homogenizationMatrix;

  this->meshIdDisp = other.meshIdDisp;
  this->dispOrder = other.dispOrder;
  this->bctype = other.bctype;

  this->leftFaces = other.leftFaces;
  this->rightFaces = other.rightFaces;
  this->topFaces = other.topFaces;
  this->bottomFaces = other.bottomFaces;
  this->frontFaces = other.frontFaces;
  this->backFaces = other.backFaces;

  this->xmin = other.xmin;
  this->xmax = other.xmax;
}

void Homogenization3DSolid::init(PointerCollection &pointers,
                                 ParameterList &parameters) {

  meshIdDisp = parameters.getIndexVal("meshIdDisp");
  dispOrder = parameters.getIndexVal("dispOrder");
  bctype = parameters.getIndexVal("bctype");

  
  auto &Logger = pointers.getSPDLogger();
  Logger.info("Initializing Homogenization3DSolid.....");
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

  this->leftFaces = getFaceNumbers(pointers, normal, this->xmin);
  this->rightFaces = getFaceNumbers(pointers, normal, this->xmax);

  Types::Vector3<prec> normal1 = {0.0, 1.0, 0.0};

  this->backFaces = getFaceNumbers(pointers, normal1, this->xmin);
  this->frontFaces = getFaceNumbers(pointers, normal1, this->xmax);

  Types::Vector3<prec> normal2 = {0.0, 0.0, 1.0};

  this->bottomFaces = getFaceNumbers(pointers, normal2, this->xmin);
  this->topFaces = getFaceNumbers(pointers, normal2, this->xmax);

  this->setDisplacementBoundaryConditions(pointers);
}

void Homogenization3DSolid::setDisplacementBoundaryConditions(
    PointerCollection &pointers) {

  Types::Vector3<indexType> dofs = {1, 1, 1};
  for (auto i : leftFaces) {
    auto face = pointers.getGeometryData()->getFace(i);
    face->setBoundaryCondition(pointers, meshIdDisp, dispOrder,
                               Geometry::ShapeFunctionTypes::H1, dofs, true);
  }
  for (auto i : rightFaces) {
    auto face = pointers.getGeometryData()->getFace(i);
    face->setBoundaryCondition(pointers, meshIdDisp, dispOrder,
                               Geometry::ShapeFunctionTypes::H1, dofs, true);
  }
  for (auto i : bottomFaces) {
    auto face = pointers.getGeometryData()->getFace(i);
    face->setBoundaryCondition(pointers, meshIdDisp, dispOrder,
                               Geometry::ShapeFunctionTypes::H1, dofs, true);
  }
  for (auto i : topFaces) {
    auto face = pointers.getGeometryData()->getFace(i);
    face->setBoundaryCondition(pointers, meshIdDisp, dispOrder,
                               Geometry::ShapeFunctionTypes::H1, dofs, true);
  }
  for (auto i : backFaces) {
    auto face = pointers.getGeometryData()->getFace(i);
    face->setBoundaryCondition(pointers, meshIdDisp, dispOrder,
                               Geometry::ShapeFunctionTypes::H1, dofs, true);
  }
  for (auto i : frontFaces) {
    auto face = pointers.getGeometryData()->getFace(i);
    face->setBoundaryCondition(pointers, meshIdDisp, dispOrder,
                               Geometry::ShapeFunctionTypes::H1, dofs, true);
  }
}

void Homogenization3DSolid::computeAMatrix(PointerCollection &pointers) {
  indexType totalEqs =
      pointers.getEquationHandler()->getNumberOfTotalEquations();
  indexType activeEqs =
      pointers.getEquationHandler()->getNumberOfActiveEquations();
  indexType inactiveIds = totalEqs - activeEqs;
  this->homogenizationMatrix.resize(inactiveIds, 6);
  this->homogenizationMatrix.setZero();
  std::vector<std::vector<indexType>> ffaces = {
      leftFaces, rightFaces, bottomFaces, topFaces, backFaces, frontFaces};
  if (bctype == 0) {
    for (auto j = 0; j < ffaces.size(); j++) {
      for (auto fNum : ffaces[j]) {
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
            this->homogenizationMatrix(posB, 1) = coor(1); // y
            this->homogenizationMatrix(posC, 2) = coor(2); // z

            this->homogenizationMatrix(posA, 3) = coor(1) * prec(0.5); // 0.5*y
            this->homogenizationMatrix(posA, 4) = coor(2) * prec(0.5); // 0.5*z
            this->homogenizationMatrix(posB, 3) = coor(0) * prec(0.5); // 0.5*x
            this->homogenizationMatrix(posB, 5) = coor(2) * prec(0.5); // 0.5*z
            this->homogenizationMatrix(posC, 4) = coor(0) * prec(0.5); // 0.5*x
            this->homogenizationMatrix(posC, 5) = coor(1) * prec(0.5); // 0.5*y
          }
        }
      }
    }
  }

  pointers.getSPDLogger().debug("Computed A-Matrix:\n{}",
                                this->homogenizationMatrix);
}

auto Homogenization3DSolid::getNumberOfStrains() -> indexType {
  return indexType(6);
}

auto Homogenization3DSolid::getAMatrix() -> Types::MatrixXX<prec> & {
  return this->homogenizationMatrix;
}

auto Homogenization3DSolid::getDv() -> prec {
  return (xmax(0) - xmin(0)) * (xmax(1) - xmin(1)) * (xmax(2) - xmin(2));
}

auto Homogenization3DSolid::getDisplacementIncrement(
    Types::VectorX<prec> &strainIncrement) -> Types::VectorX<prec> {
  Types::VectorX<prec> disp = this->homogenizationMatrix * strainIncrement;
  return disp;
}

void Homogenization3DSolid::toFile(PointerCollection &pointers,
                                   std::ofstream &out) {
  HomogenizationBase::toFile(pointers, out);
  writeScalar(out, meshIdDisp);
  writeScalar(out, dispOrder);
  writeScalar(out, bctype);
  writeStdVector(out, leftFaces);
  writeStdVector(out, rightFaces);
  writeStdVector(out, bottomFaces);
  writeStdVector(out, topFaces);
  writeStdVector(out, backFaces);
  writeStdVector(out, frontFaces);

  writeEigenMatrix(out, homogenizationMatrix);
  writeEigenMatrix(out, xmin);
  writeEigenMatrix(out, xmax);
}

void Homogenization3DSolid::fromFile(PointerCollection &pointers,
                                     std::ifstream &in) {
  HomogenizationBase::fromFile(pointers, in);
  readScalar(in, meshIdDisp);
  readScalar(in, dispOrder);
  readScalar(in, bctype);
  readStdVector(in, leftFaces);
  readStdVector(in, rightFaces);
  readStdVector(in, bottomFaces);
  readStdVector(in, topFaces);
  readStdVector(in, backFaces);
  readStdVector(in, frontFaces);

  readEigenMatrix(in, homogenizationMatrix);
  readEigenMatrix(in, xmin);
  readEigenMatrix(in, xmax);
}

} // namespace HierAMuS