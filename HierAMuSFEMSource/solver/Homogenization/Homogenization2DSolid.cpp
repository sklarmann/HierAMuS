// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "Homogenization2DSolid.h"
#include "datatypes.h"
#include "solver/Homogenization/Homogenization2DSolid.h"


// Equations 
#include "EquationHandler.h"


#include "geometry/GeometryData.h"
#include "geometry/Edges/EdgesData.h"
#include "geometry/VertexData.h"

#include "control/BinaryWrite.h"

#include "spdlog/fmt/ostr.h"

namespace HierAMuS {
Homogenization2DSolid::Homogenization2DSolid()
    : xmin(0.0), xmax(0.0), ymin(0.0), ymax(0.0) {}

Homogenization2DSolid::Homogenization2DSolid(
    const Homogenization2DSolid &other) : HomogenizationBase(other)
{
  this->homogenizationMatrix = other.homogenizationMatrix;
  this->meshIdDisp = other.meshIdDisp;
  this->dispOrder = other.dispOrder;
  this->bctype = other.bctype;
  this->leftedges = other.leftedges;
  this->rightedges = other.rightedges;
  this->bottomedges = other.bottomedges;
  this->topedges = other.topedges;

  this->xmin = other.xmin;
  this->xmax = other.xmax;
  this->ymin = other.ymin;
  this->ymax = other.ymax;
}

void Homogenization2DSolid::init(PointerCollection &pointers,
                                 ParameterList &parameters) {
  
  meshIdDisp = parameters.getIndexVal("meshIdDisp");
  bctype = parameters.getIndexVal("bctype");
  dispOrder = parameters.getIndexVal("dispOrder");

  auto &Logger = pointers.getSPDLogger();
  Logger.info("Initializing Homogenization2DSolid.....");
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

  auto getEdgeNumbers = [](PointerCollection &pointers, Types::Vector3<prec> &normal, Types::Vector3<prec> & point)
  {
    auto edges =
        pointers.getGeometryData()->getEdgesInPlane(normal, point);
    Types::MatrixXX<indexType> edgeNums;
    indexType nEdges = edges.size();
    edgeNums.resize(nEdges, 1);
    edgeNums.setZero();
    for (indexType i=0;i<nEdges;++i)
    {
      edgeNums(i, 0) = edges[i]->getId();
    }
    return edgeNums;
  };

  // getEdges
  Types::Vector3<prec> normal, point;
  normal.setZero();
  point.setZero();
  Types::Vector3<prec> xMin = pointers.getGeometryData()->getxMin();
  Types::Vector3<prec> xMax = pointers.getGeometryData()->getxMax();
  this->xmin = xMin(0);
  this->xmax = xMax(0);
  this->ymin = xMin(1);
  this->ymax = xMax(1);

  normal(0) = prec(1);
  std::cout << leftedges << std::endl;
  leftedges = getEdgeNumbers(pointers, normal, xMin);
  std::cout << "new left edges:\n" << leftedges << std::endl;
  rightedges = getEdgeNumbers(pointers, normal, xMax);
  normal(0) = prec(0);
  normal(1) = prec(1);
  bottomedges = getEdgeNumbers(pointers, normal, xMin);
  topedges = getEdgeNumbers(pointers, normal, xMax);
  


  if (bctype == 0) {
    setDisplacementBoundaryConditions(pointers, meshIdDisp, dispOrder,
                                      leftedges, rightedges, bottomedges,
                                      topedges);
  }
}

void Homogenization2DSolid::setEdgesBC(PointerCollection &pointers,
                                       indexType meshIdDisp,
                                       indexType dispOrder,
                                       Types::MatrixXX<indexType> &edges) {
  auto geoData = pointers.getGeometryData();

  Types::Vector3<indexType> bc = {1, 1, 1};
  for (auto i = 0; i < edges.cols(); ++i) {
    for (auto j = 0; j < edges.rows(); ++j) {
      indexType ednum = edges(j, i);
      auto &Edge = geoData->getEdgeData(ednum);
      Edge.setBoundaryCondition(meshIdDisp, dispOrder, Geometry::ShapeFunctionTypes::H1,
                                bc, true);
    }
  }
}

void Homogenization2DSolid::setDisplacementBoundaryConditions(
    PointerCollection &pointers, indexType meshIdDisp, indexType dispOrder,
    Types::MatrixXX<indexType> &leftedges,
    Types::MatrixXX<indexType> &rightedges,
    Types::MatrixXX<indexType> &bottomedges,
    Types::MatrixXX<indexType> &topedges) {

  setEdgesBC(pointers, meshIdDisp, dispOrder, leftedges);
  setEdgesBC(pointers, meshIdDisp, dispOrder, rightedges);
  setEdgesBC(pointers, meshIdDisp, dispOrder, topedges);
  setEdgesBC(pointers, meshIdDisp, dispOrder, bottomedges);
}

void Homogenization2DSolid::computeAMatrix(PointerCollection &pointers) {
  indexType totalEqs = pointers.getEquationHandler()->getNumberOfTotalEquations();
  indexType activeEqs = pointers.getEquationHandler()->getNumberOfActiveEquations();
  indexType inactiveIds = totalEqs - activeEqs;

  this->homogenizationMatrix.resize(inactiveIds, 3);
  this->homogenizationMatrix.setZero();

  edgeEntriesAMatrix(pointers, meshIdDisp, dispOrder, leftedges);
  edgeEntriesAMatrix(pointers, meshIdDisp, dispOrder, rightedges);
  edgeEntriesAMatrix(pointers, meshIdDisp, dispOrder, bottomedges);
  edgeEntriesAMatrix(pointers, meshIdDisp, dispOrder, topedges);


  pointers.getSPDLogger().trace("Computed A-Matrix:\n{}",this->homogenizationMatrix);
}

auto Homogenization2DSolid::getNumberOfStrains() -> indexType {
  return indexType(3);
}

auto Homogenization2DSolid::getAMatrix() -> Types::MatrixXX<prec> & {
  return this->homogenizationMatrix;
}

auto Homogenization2DSolid::getDv() -> prec {
  prec dlx = xmax - xmin;
  prec dly = ymax - ymin;
  return dlx * dly;
}

void Homogenization2DSolid::edgeEntriesAMatrix(
    PointerCollection &pointers, indexType meshIdDisp, indexType dispOrder,
    Types::MatrixXX<indexType> &edges) {
  auto geoData = pointers.getGeometryData();
  for (auto i = 0; i < edges.rows(); ++i) {
    for (auto j = 0; j < edges.cols(); ++j) {
      indexType ednum = edges(i, j);
      auto &Edge = geoData->getEdgeData(ednum);
      for (auto k = 0; k < Edge.getNumberOfVerts(); ++k) {
        auto &Vert = *Edge.getVertex(k);
        auto Nodes = Vert.getNodesOfSet(meshIdDisp);
        auto coor = Vert.getCoordinates();
        auto Dofs = Nodes[0]->getDegreesOfFreedom();

        indexType posx = Dofs[0]->getEqId();
        indexType posy = Dofs[1]->getEqId();
        this->homogenizationMatrix(posx, 0) = coor(0);
        this->homogenizationMatrix(posx, 2) = coor(1) * prec(0.5);

        this->homogenizationMatrix(posy, 1) = coor(1);
        this->homogenizationMatrix(posy, 2) = coor(0) * prec(0.5);
        if (coor(0) > xmax)
          xmax = coor(0);
        if (coor(0) < xmin)
          xmin = coor(0);
        if (coor(1) > ymax)
          ymax = coor(1);
        if (coor(1) < ymin)
          ymin = coor(1);
      }
    }
  }
}

auto Homogenization2DSolid::getDisplacementIncrement(
    Types::VectorX<prec> &strainIncrement) -> Types::VectorX<prec> {
  Types::VectorX<prec> disp =
      this->homogenizationMatrix * strainIncrement;
  return disp;
}

void Homogenization2DSolid::toFile(PointerCollection &pointers,
                                   std::ofstream &out) {
  HomogenizationBase::toFile(pointers, out);
  writeScalar(out, meshIdDisp);
  writeScalar(out, dispOrder);
  writeScalar(out, bctype);
  writeScalar(out, xmin);
  writeScalar(out, xmax);
  writeScalar(out, ymin);
  writeScalar(out, ymax);

  writeEigenMatrix(out, homogenizationMatrix);
  writeEigenMatrix(out, leftedges);
  writeEigenMatrix(out, rightedges);
  writeEigenMatrix(out, bottomedges);
  writeEigenMatrix(out, topedges);
}

void Homogenization2DSolid::fromFile(PointerCollection &pointers,
                                     std::ifstream &in) {
  HomogenizationBase::fromFile(pointers, in);
  readScalar(in, meshIdDisp);
  readScalar(in, dispOrder);
  readScalar(in, bctype);
  readScalar(in, xmin);
  readScalar(in, xmax);
  readScalar(in, ymin);
  readScalar(in, ymax);

  readEigenMatrix(in, homogenizationMatrix);
  readEigenMatrix(in, leftedges);
  readEigenMatrix(in, rightedges);
  readEigenMatrix(in, bottomedges);
  readEigenMatrix(in, topedges);
}

} // namespace HierAMuS