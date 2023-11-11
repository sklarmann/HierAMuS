// Copyright 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "Homogenization3DThermoMechBeam.h"
#include "MatrixTypes.h"
#include "datatypes.h"
#include "solver/Constraints/GeneralLink.h"
#include "solver/Homogenization/Homogenization3DThermoMechBeam.h"
#include "solver/GenericSolutionState.h"

//Equations
#include "EquationHandler.h"

#include "geometry/GeometryData.h"
#include "geometry/GeometryTypes.h"
#include "geometry/VertexData.h"
#include "geometry/Faces/FacesData.h"

#include "control/BinaryWrite.h"

#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/bundled/format.h"

#include <algorithm>

namespace HierAMuS {
Homogenization3DThermoMechBeam::Homogenization3DThermoMechBeam()
    : xmin({0.0, 0.0, 0.0}), xmax({0.0, 0.0, 0.0}) {}

Homogenization3DThermoMechBeam::Homogenization3DThermoMechBeam(
    const Homogenization3DThermoMechBeam &other)
    : HomogenizationBase(other) {
  this->homogenizationMatrix = other.homogenizationMatrix;
  this->homogenizationMatrixDisp = other.homogenizationMatrixDisp;
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

  auto &Logger = pointers.getSPDLogger();

  if (bctype == 0) {
    Logger.info("\n{:-<100}\n"
                "Initializing Homogenization3DThermoMechBeam.....\n"
                "   Mesh Id for displacements:                 {:>12}\n"
                "   Order of approximation of displacements:   {:>12}\n"
                "   Mesh Id for temperature:                   {:>12}\n"
                "   Order of approximation of temperature:     {:>12}\n"
                "   Using displacement boundary conditions."
                "\n{:-<100}\n",
                "", meshIdDisp, dispOrder, meshIdTemp, tempOrder, "");

  } else if (bctype == 1) {
    Logger.info(
        "\n{:-<100}\n"
        "Initializing Homogenization3DThermoMechBeam.....\n"
        "   Mesh Id for displacements:                 {:>12}\n"
        "   Order of approximation of displacements:   {:>12}\n"
        "   Mesh Id for temperature:                   {:>12}\n"
        "   Order of approximation of temperature:     {:>12}\n"
        "   Using temperature and periodic displacement boundary conditions."
        "\n{:-<100}\n",
        "", meshIdDisp, dispOrder, meshIdTemp, tempOrder, "");
  } else if (bctype == 2) {
    Logger.info("\n{:-<100}\n"
                "Initializing Homogenization3DThermoMechBeam.....\n"
                "   Mesh Id for displacements:                 {:>12}\n"
                "   Order of approximation of displacements:   {:>12}\n"
                "   Mesh Id for temperature:                   {:>12}\n"
                "   Order of approximation of temperature:     {:>12}\n"
                "   Using periodic temperature and periodic displacement "
                "boundary conditions."
                "\n{:-<100}\n",
                "", meshIdDisp, dispOrder, meshIdTemp, tempOrder, "");
  } else {

    Logger.warn("\n{:-<100}\n"
                "Initializing Homogenization3DThermoMechBeam.....\n"
                "   Mesh Id for displacements:                 {:>12}\n"
                "   Order of approximation of displacements:   {:>12}\n"
                "   Mesh Id for temperature:                   {:>12}\n"
                "   Order of approximation of temperature:     {:>12}\n"
                "   Specified unsupported boundary condition type, using "
                "displacement boundary conditions."
                "\n{:-<100}\n",
                "", meshIdDisp, dispOrder, meshIdTemp, tempOrder, "");
    bctype = 0;
  }

  this->computeGeometryParameters(pointers);

  if (bctype == 0) {
    this->setDisplacementBoundaryConditions(pointers);
  } else if (bctype == 3) {
    this->setDisplacementBoundaryConditions2(pointers);
  } else if (bctype == 1 || bctype == 2) {
    pointers.getGeometryData()->sortReorientFacesPeriodicBC(
        this->leftFacesMaster, this->rightFacesSlave);
    pointers.getGeometryData()->sortReorientFacesPeriodicBC(
        this->bottomFacesMaster, this->topFacesSlave);
    pointers.getGeometryData()->sortReorientFacesPeriodicBC(
        this->backFacesMaster, this->frontFacesSlave);
    Logger.debug(
        "Arranged master and slave faces for periodic boundary conditions\n"
        "Master faces: {}\n"
        "Slave faces: {}\n",
        fmt::join(leftFacesMaster, " "), fmt::join(rightFacesSlave, " "));
    this->setPeriodicBoundaryConditions(pointers);
  }
}

void Homogenization3DThermoMechBeam::computeGeometryParameters(
    PointerCollection &pointers) {

  // Setting up the faces
  auto getFaceNumbers = [](PointerCollection &pointers,
                           Types::Vector3<prec> &normal,
                           Types::Vector3<prec> &point) {
    auto faces =
        pointers.getGeometryData()->getFacesInPlane(normal, point);
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

  normal = {0.0, 1.0, 0.0};
  this->backFacesMaster = getFaceNumbers(pointers, normal, this->xmin);
  this->frontFacesSlave = getFaceNumbers(pointers, normal, this->xmax);

  normal = {0.0, 0.0, 1.0};
  this->bottomFacesMaster = getFaceNumbers(pointers, normal, this->xmin);
  this->topFacesSlave = getFaceNumbers(pointers, normal, this->xmax);

  this->v1 =
      pointers.getGeometryData()
          ->getVertexClosestTo({this->xmin[0], this->xmin[1], this->xmin[2]})
          .getId();
  this->v2 =
      pointers.getGeometryData()
          ->getVertexClosestTo({this->xmax[0], this->xmin[1], this->xmin[2]})
          .getId();
  this->v3 =
      pointers.getGeometryData()
          ->getVertexClosestTo({this->xmax[0], this->xmax[1], this->xmin[2]})
          .getId();
  this->v4 =
      pointers.getGeometryData()
          ->getVertexClosestTo({this->xmin[0], this->xmax[1], this->xmin[2]})
          .getId();
  this->v5 =
      pointers.getGeometryData()
          ->getVertexClosestTo({this->xmin[0], this->xmin[1], this->xmax[2]})
          .getId();
  this->v6 =
      pointers.getGeometryData()
          ->getVertexClosestTo({this->xmax[0], this->xmin[1], this->xmax[2]})
          .getId();
  this->v7 =
      pointers.getGeometryData()
          ->getVertexClosestTo({this->xmax[0], this->xmax[1], this->xmax[2]})
          .getId();
  this->v8 =
      pointers.getGeometryData()
          ->getVertexClosestTo({this->xmin[0], this->xmax[1], this->xmax[2]})
          .getId();

  auto getEdgesDir = [](PointerCollection &pointers,
                        const Types::Vector3<prec> &normal,
                        indexType startvertex, indexType endvertex) {
    indexType csv = startvertex;
    std::vector<indexType> edges;
    while (csv != endvertex) {
      auto &V = pointers.getGeometryData()->getVertexData(csv);
      auto enums = V.getConnectedEdges();
      for (auto i : enums) {
        auto &edge = pointers.getGeometryData()->getEdgeData(i);
        IntegrationPoint ip;
        ip.xi = 0.0;
        Types::Vector3<prec> A1 = edge.getA1Vector(ip);
        if (abs(A1.dot(normal)) >= prec(0.999)) {
          if (std::find(edges.begin(), edges.end(), i) == edges.end()) {
            edges.push_back(i);
            auto eV1 = edge.getVertexNumber(0);
            auto eV2 = edge.getVertexNumber(1);
            if (eV1 == csv) {
              csv = eV2;
            } else {
              csv = eV1;
              edge.flip();
            }
          }
        }
      }
    }
    return edges;
  };

  xEdges1 = getEdgesDir(pointers, {1, 0, 0}, v1, v2);
  xEdges2 = getEdgesDir(pointers, {1, 0, 0}, v5, v6);
  xEdges3 = getEdgesDir(pointers, {1, 0, 0}, v8, v7);
  xEdges4 = getEdgesDir(pointers, {1, 0, 0}, v4, v3);

  yEdges1 = getEdgesDir(pointers, {0, 1, 0}, v1, v4);
  yEdges2 = getEdgesDir(pointers, {0, 1, 0}, v2, v3);
  yEdges3 = getEdgesDir(pointers, {0, 1, 0}, v5, v8);
  yEdges4 = getEdgesDir(pointers, {0, 1, 0}, v6, v7);

  zEdges1 = getEdgesDir(pointers, {0, 0, 1}, v1, v5);
  zEdges2 = getEdgesDir(pointers, {0, 0, 1}, v2, v6);
  zEdges3 = getEdgesDir(pointers, {0, 0, 1}, v3, v7);
  zEdges4 = getEdgesDir(pointers, {0, 0, 1}, v4, v8);

  auto &Logger = pointers.getSPDLogger();

  Logger.debug(
      "Computed geometry information for Homogenization3DThermoMechBeam...\n"
      "  With left faces (Master):   {}\n"
      "  With right faces (Slave):   {}\n"
      "  With back faces (Master):   {}\n"
      "  With fron faces (Slave):    {}\n"
      "  With bottom faces (Master): {}\n"
      "  With top faces (Slave):     {}\n"
      "  x-Edges1 (Master):          {}\n"
      "  x-Edges2 (Slave):           {}\n"
      "  x-Edges3 (Slave):           {}\n"
      "  x-Edges4 (Slave):           {}\n"
      "  y-Edges1 (Master):          {}\n"
      "  y-Edges2 (Slave):           {}\n"
      "  y-Edges3 (Slave):           {}\n"
      "  y-Edges4 (Slave):           {}\n"
      "  z-Edges1 (Master):          {}\n"
      "  z-Edges2 (Slave):           {}\n"
      "  z-Edges3 (Slave):           {}\n"
      "  z-Edges4 (Slave):           {}\n"
      "  Corner vertices:            {}\n"
      "  ",
      fmt::join(this->leftFacesMaster, ", "),
      fmt::join(this->rightFacesSlave, ", "),
      fmt::join(this->backFacesMaster, ", "),
      fmt::join(this->frontFacesSlave, ", "),
      fmt::join(this->bottomFacesMaster, ", "),
      fmt::join(this->topFacesSlave, ", "), fmt::join(this->xEdges1, ", "),
      fmt::join(this->xEdges2, ", "), fmt::join(this->xEdges3, ", "),
      fmt::join(this->xEdges4, ", "), fmt::join(this->yEdges1, ", "),
      fmt::join(this->yEdges2, ", "), fmt::join(this->yEdges3, ", "),
      fmt::join(this->yEdges4, ", "), fmt::join(this->zEdges1, ", "),
      fmt::join(this->zEdges2, ", "), fmt::join(this->zEdges3, ", "),
      fmt::join(this->zEdges4, ", "),
      fmt::join(std::vector<indexType>({v1, v2, v3, v4, v5, v6, v7, v8}), ", ")

  );
}

void Homogenization3DThermoMechBeam::setDisplacementBoundaryConditions(
    PointerCollection &pointers) {

  Types::Vector3<indexType> dofs = {1, 1, 1};
  for (const auto &v :
       {std::cref(leftFacesMaster), std::cref(rightFacesSlave)}) {
    for (auto i : v.get()) {
      auto face = pointers.getGeometryData()->getFaceData(i);
      face->setBoundaryCondition(meshIdDisp, dispOrder,
                                 Geometry::ShapeFunctionTypes::H1, dofs, true);
      face->setBoundaryCondition(meshIdTemp, tempOrder,
                                 Geometry::ShapeFunctionTypes::H1, dofs, true);
    }
  }

  for (const auto &v :
       {std::cref(bottomFacesMaster), std::cref(topFacesSlave),
        std::cref(backFacesMaster), std::cref(frontFacesSlave)}) {
    for (auto i : v.get()) {
      auto face = pointers.getGeometryData()->getFaceData(i);
      face->setBoundaryCondition(meshIdTemp, tempOrder,
                                 Geometry::ShapeFunctionTypes::H1, dofs, true);
    }
  }
}

void Homogenization3DThermoMechBeam::setDisplacementBoundaryConditions2(
    PointerCollection &pointers) {
  Types::Vector3<indexType> dofs = {1, 1, 1};
  for (const auto &v :
       {std::cref(leftFacesMaster), std::cref(rightFacesSlave)}) {
    for (auto i : v.get()) {
      auto face = pointers.getGeometryData()->getFaceData(i);
      face->setBoundaryCondition(meshIdDisp, dispOrder,
                                 Geometry::ShapeFunctionTypes::H1, dofs, true);
      face->setBoundaryCondition(meshIdTemp, tempOrder,
                                 Geometry::ShapeFunctionTypes::H1, dofs, true);
    }
  }
}

void Homogenization3DThermoMechBeam::setPeriodicBoundaryConditions(
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

    // Types::Vector3<indexType> dofs = {1, 1, 1};
    // for (const auto &v :
    //      {std::cref(bottomFacesMaster), std::cref(topFacesSlave),
    //       std::cref(backFacesMaster), std::cref(frontFacesSlave),
    //       std::cref(leftFacesMaster), std::cref(rightFacesSlave)}) {
    //   for (auto i : v.get()) {
    //     auto face = pointers.getGeometryData()->getFace(i);
    //     face->setBoundaryCondition(pointers, meshIdTemp, tempOrder,
    //                                Geometry::ShapeFunctionTypes::H1, dofs,
    //                                true);
    //   }
    // }
    Types::Vector3<indexType> dofs = {1, 1, 1};
    for (const auto &v :
         {std::cref(leftFacesMaster), std::cref(rightFacesSlave),
          std::cref(topFacesSlave), std::cref(bottomFacesMaster),
          std::cref(backFacesMaster), std::cref(frontFacesSlave)}) {
      for (auto i : v.get()) {
        auto face = pointers.getGeometryData()->getFaceData(i);
        face->setBoundaryCondition(meshIdTemp, tempOrder,
                                   Geometry::ShapeFunctionTypes::H1, dofs,
                                   true);
      }
    }
  } else if (bctype == 2) {
    // Displacements
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

    // Thermal vertices
    for (auto i : {v1, v2, v3, v4, v5, v6, v7, v8}) {
      auto &vert = pointers.getGeometryData()->getVertexData(i);
      Types::Vector3<indexType> bc = {1, 1, 1};
      vert.setBoundaryCondition(meshIdTemp, tempOrder,
                                Geometry::ShapeFunctionTypes::H1, bc, true);
    }
    // Thermal edges x direction
    for (indexType i = 0; i < 1; ++i) {
      pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
          pointers, Geometry::GeometryTypes::Edges, xEdges1, xEdges2,
          meshIdTemp, tempOrder, i, i, 1, 0, true);
      pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
          pointers, Geometry::GeometryTypes::Edges, xEdges1, xEdges3,
          meshIdTemp, tempOrder, i, i, 1, 0, true);
      pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
          pointers, Geometry::GeometryTypes::Edges, xEdges1, xEdges4,
          meshIdTemp, tempOrder, i, i, 1, 0, true);
    }
    // Thermal edges y direction
    for (indexType i = 0; i < 1; ++i) {
      pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
          pointers, Geometry::GeometryTypes::Edges, yEdges1, yEdges2,
          meshIdTemp, tempOrder, i, i, 1, 0, true);
      pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
          pointers, Geometry::GeometryTypes::Edges, yEdges1, yEdges3,
          meshIdTemp, tempOrder, i, i, 1, 0, true);
      pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
          pointers, Geometry::GeometryTypes::Edges, yEdges1, yEdges4,
          meshIdTemp, tempOrder, i, i, 1, 0, true);
    }
    // Thermal edges z direction
    for (indexType i = 0; i < 1; ++i) {
      pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
          pointers, Geometry::GeometryTypes::Edges, zEdges1, zEdges2,
          meshIdTemp, tempOrder, i, i, 1, 0, true);
      pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
          pointers, Geometry::GeometryTypes::Edges, zEdges1, zEdges3,
          meshIdTemp, tempOrder, i, i, 1, 0, true);
      pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
          pointers, Geometry::GeometryTypes::Edges, zEdges1, zEdges4,
          meshIdTemp, tempOrder, i, i, 1, 0, true);
    }

    // Faces
    for (indexType i = 0; i < 3; ++i) {
      pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
          pointers, Geometry::GeometryTypes::Faces, this->leftFacesMaster,
          this->rightFacesSlave, this->meshIdTemp, this->tempOrder, i, i, 1, 0,
          true);
      pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
          pointers, Geometry::GeometryTypes::Faces, this->bottomFacesMaster,
          this->topFacesSlave, this->meshIdTemp, this->tempOrder, i, i, 1, 0,
          true);
      pointers.getSolutionState()->getConstraintHandler().GeneralLinkGeo(
          pointers, Geometry::GeometryTypes::Faces, this->backFacesMaster,
          this->frontFacesSlave, this->meshIdTemp, this->tempOrder, i, i, 1, 0,
          true);
    }
  }
}

void Homogenization3DThermoMechBeam::computeAMatrix(
    PointerCollection &pointers) {
  indexType totalEqs =
      pointers.getEquationHandler()->getNumberOfTotalEquations();
  indexType activeEqs =
      pointers.getEquationHandler()->getNumberOfActiveEquations();
  indexType inactiveIds = totalEqs - activeEqs;

  this->homogenizationMatrix.resize(inactiveIds, 12);
  this->homogenizationMatrix.setZero();

  if (bctype == 0) {
    for (const auto &ff :
         {std::cref(leftFacesMaster), std::cref(rightFacesSlave)}) {
      for (auto fNum : ff.get()) {
        auto face = pointers.getGeometryData()->getFaceData(fNum);
        indexType numV = face->getNumberOfVerts();
        for (indexType i = 0; i < numV; ++i) {
          auto vert = face->getVertex(i);
          auto Nodes = vert->getNodesOfSet(meshIdDisp);
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
          auto NodesT = vert->getNodesOfSet(meshIdTemp);
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
         {std::cref(bottomFacesMaster), std::cref(topFacesSlave),
          std::cref(frontFacesSlave), std::cref(backFacesMaster)}) {
      for (auto fNum : ff.get()) {
        auto face = pointers.getGeometryData()->getFaceData(fNum);
        indexType numV = face->getNumberOfVerts();
        for (indexType i = 0; i < numV; ++i) {
          auto vert = face->getVertex(i);
          auto Nodes = vert->getNodesOfSet(meshIdTemp);
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
        auto face = pointers.getGeometryData()->getFaceData(fNum);
        indexType numV = face->getNumberOfVerts();
        for (indexType i = 0; i < numV; ++i) {
          auto vert = face->getVertex(i);
          auto Nodes = vert->getNodesOfSet(meshIdDisp);
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
          auto NodesT = vert->getNodesOfSet(meshIdTemp);
          for (auto node : NodesT) {
            indexType posA = node->getDegreeOfFreedom(0).getEqId();
            this->homogenizationMatrix(posA, 6) = prec(1.0);           // x
            this->homogenizationMatrix(posA, 7) = coor(2);             // x*z
            this->homogenizationMatrix(posA, 8) = -coor(1);            // -x*y
            this->homogenizationMatrix(posA, 9) = coor(0);             // x
            this->homogenizationMatrix(posA, 10) = coor(0) * coor(2);  // x*z
            this->homogenizationMatrix(posA, 11) = -coor(0) * coor(1); // -x*y
          }
        }
      }
    }

  } else if (bctype == 1) {
    for (indexType i = 0; i < static_cast<indexType>(leftFacesMaster.size()); ++i) {
      auto mFace = pointers.getGeometryData()->getFaceData(leftFacesMaster[i]);
      auto sFace = pointers.getGeometryData()->getFaceData(rightFacesSlave[i]);
      indexType numv = mFace->getNumberOfVerts();
      for (indexType lv = 0; lv < numv; ++lv) {
        auto mVert = mFace->getVertex(lv);
        auto sVert = sFace->getVertex(lv);
        Types::Vector3<prec> dx =
            sVert->getCoordinates() - mVert->getCoordinates();

        Types::Vector3<prec> coor = sVert->getCoordinates();
        auto sNodes = sVert->getNodesOfSet(meshIdDisp);
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
    for (const auto &ff :
         {std::cref(leftFacesMaster), std::cref(rightFacesSlave)}) {
      for (auto fNum : ff.get()) {
        auto face = pointers.getGeometryData()->getFaceData(fNum);
        indexType numV = face->getNumberOfVerts();
        for (indexType i = 0; i < numV; ++i) {
          auto vert = face->getVertex(i);
          auto Nodes = vert->getNodesOfSet(meshIdTemp);
          auto coor = vert->getCoordinates();
          for (auto node : Nodes) {
            indexType posA = node->getDegreeOfFreedom(0).getEqId();
            this->homogenizationMatrix(posA, 6) = prec(1.0);           // x
            this->homogenizationMatrix(posA, 7) = coor(2);//*prec(2); // x*z
            this->homogenizationMatrix(posA, 8) = -coor(1);//*prec(2); // -x*y
            this->homogenizationMatrix(posA, 9) = coor(0);             // x
            this->homogenizationMatrix(posA, 10) = coor(0) * coor(2);  // x*z
            this->homogenizationMatrix(posA, 11) = -coor(0) * coor(1); // -x*y
          }
        }
      }
    }
    for (const auto &ff :
         {std::cref(topFacesSlave), std::cref(bottomFacesMaster),
          std::cref(frontFacesSlave), std::cref(backFacesMaster)}) {
      for (auto fNum : ff.get()) {
        auto face = pointers.getGeometryData()->getFaceData(fNum);
        indexType numV = face->getNumberOfVerts();
        for (indexType i = 0; i < numV; ++i) {
          auto vert = face->getVertex(i);
          auto Nodes = vert->getNodesOfSet(meshIdTemp);
          auto coor = vert->getCoordinates();
          for (auto node : Nodes) {
            indexType posA = node->getDegreeOfFreedom(0).getEqId();
            this->homogenizationMatrix(posA, 6) = prec(1.0);           // x
            this->homogenizationMatrix(posA, 7) = coor(2);             // x*z
            this->homogenizationMatrix(posA, 8) = -coor(1);            // -x*y
            this->homogenizationMatrix(posA, 9) = coor(0);             // x
            this->homogenizationMatrix(posA, 10) = coor(0) * coor(2);  // x*z
            this->homogenizationMatrix(posA, 11) = -coor(0) * coor(1); // -x*y
          }
        }
      }
    }
  } else if (bctype == 2) {
    this->homogenizationMatrixDisp.resize(inactiveIds, 12);
    this->homogenizationMatrixDisp.setZero();

    // Periodic displacements
    for (indexType i = 0; i < static_cast<indexType>(leftFacesMaster.size()); ++i) {
      auto mFace = pointers.getGeometryData()->getFaceData(leftFacesMaster[i]);
      auto sFace = pointers.getGeometryData()->getFaceData(rightFacesSlave[i]);
      indexType numv = mFace->getNumberOfVerts();
      for (indexType lv = 0; lv < numv; ++lv) {
        auto mVert = mFace->getVertex(lv);
        auto sVert = sFace->getVertex(lv);
        Types::Vector3<prec> dx =
            sVert->getCoordinates() - mVert->getCoordinates();

        Types::Vector3<prec> coor = sVert->getCoordinates();
        auto sNodes = sVert->getNodesOfSet(meshIdDisp);
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
    // Thermal periodic  left right faces
    for (indexType i = 0; i < static_cast<indexType>(leftFacesMaster.size()); ++i) {
      auto mFace = pointers.getGeometryData()->getFaceData(leftFacesMaster[i]);
      auto sFace = pointers.getGeometryData()->getFaceData(rightFacesSlave[i]);
      indexType numv = mFace->getNumberOfVerts();
      for (indexType lv = 0; lv < numv; ++lv) {
        auto mVert = mFace->getVertex(lv);
        auto sVert = sFace->getVertex(lv);
        Types::Vector3<prec> dx =
            sVert->getCoordinates() - mVert->getCoordinates();

        Types::Vector3<prec> coor = sVert->getCoordinates();
        auto sNodes = sVert->getNodesOfSet(meshIdTemp);
        indexType posA = sNodes[0]->getDegreeOfFreedom(0).getEqId();

        this->homogenizationMatrix(posA, 7) = dx(2);             // x
        this->homogenizationMatrix(posA, 8) = -dx(1);            // x*z
        this->homogenizationMatrix(posA, 9) = dx(0);             // -x*y
        this->homogenizationMatrix(posA, 10) = dx(0) * coor(2);  // x*z
        this->homogenizationMatrix(posA, 11) = -dx(0) * coor(1); // -x*y
      }
    }
    // Thermal periodic  top bottom faces
    for (indexType i = 0; i < static_cast<indexType>(bottomFacesMaster.size()); ++i) {
      auto mFace = pointers.getGeometryData()->getFaceData(bottomFacesMaster[i]);
      auto sFace = pointers.getGeometryData()->getFaceData(topFacesSlave[i]);
      indexType numv = mFace->getNumberOfVerts();
      for (indexType lv = 0; lv < numv; ++lv) {
        auto mVert = mFace->getVertex(lv);
        auto sVert = sFace->getVertex(lv);
        Types::Vector3<prec> dx =
            sVert->getCoordinates() - mVert->getCoordinates();

        Types::Vector3<prec> coor = sVert->getCoordinates();
        auto sNodes = sVert->getNodesOfSet(meshIdTemp);
        indexType posA = sNodes[0]->getDegreeOfFreedom(0).getEqId();

        this->homogenizationMatrix(posA, 7) = dx(2);            // x
        this->homogenizationMatrix(posA, 10) = coor(0) * dx(2); // x*z
      }
    }
    // Thermal periodic  front back faces
    for (indexType i = 0; i < static_cast<indexType>(backFacesMaster.size()); ++i) {
      auto mFace = pointers.getGeometryData()->getFaceData(backFacesMaster[i]);
      auto sFace = pointers.getGeometryData()->getFaceData(frontFacesSlave[i]);
      indexType numv = mFace->getNumberOfVerts();
      for (indexType lv = 0; lv < numv; ++lv) {
        auto mVert = mFace->getVertex(lv);
        auto sVert = sFace->getVertex(lv);
        Types::Vector3<prec> dx =
            sVert->getCoordinates() - mVert->getCoordinates();

        Types::Vector3<prec> coor = sVert->getCoordinates();
        auto sNodes = sVert->getNodesOfSet(meshIdTemp);
        indexType posA = sNodes[0]->getDegreeOfFreedom(0).getEqId();

        
        this->homogenizationMatrix(posA, 8) = -dx(1);            // x
        this->homogenizationMatrix(posA, 11) = -coor(0) * dx(1); // x*z
      }
    }

    // Thermal periodic xEdges
    for (auto &sEdges : {xEdges2, xEdges3, xEdges4}) {
      for (auto i = 0; i < xEdges1.size(); ++i) {
        auto &mEdge = pointers.getGeometryData()->getEdgeData(xEdges1[i]);
        auto &sEdge1 = pointers.getGeometryData()->getEdgeData(sEdges[i]);

        for (auto i = 0; i < 2; ++i) {
          auto &mV = *mEdge.getVertex(i);
          auto &sV1 = *sEdge1.getVertex(i);
          Types::Vector3<prec> dx = sV1.getCoordinates() - mV.getCoordinates();
          Types::Vector3<prec> coor = sV1.getCoordinates();
          auto sNodes = sV1.getNodesOfSet(meshIdTemp);
          indexType posA = sNodes[0]->getDegreeOfFreedom(0).getEqId();

          this->homogenizationMatrix(posA, 7) = dx(2);             // x
          this->homogenizationMatrix(posA, 8) = -dx(1);            // x
          this->homogenizationMatrix(posA, 9) = dx(0);             // x
          this->homogenizationMatrix(posA, 10) = coor(0) * dx(2);  // x
          this->homogenizationMatrix(posA, 11) = -coor(0) * dx(1); // x
        }
      }
    }

    // Thermal periodic yEdges
    for (auto &sEdges : {yEdges2, yEdges3, yEdges4}) {
      for (auto i = 0; i < yEdges1.size(); ++i) {
        auto &mEdge = pointers.getGeometryData()->getEdgeData(yEdges1[i]);
        auto &sEdge1 = pointers.getGeometryData()->getEdgeData(sEdges[i]);

        for (auto i = 0; i < 2; ++i) {
          auto &mV = *mEdge.getVertex(i);
          auto &sV1 = *sEdge1.getVertex(i);
          Types::Vector3<prec> dx = sV1.getCoordinates() - mV.getCoordinates();
          Types::Vector3<prec> coor = sV1.getCoordinates();
          auto sNodes = sV1.getNodesOfSet(meshIdTemp);
          indexType posA = sNodes[0]->getDegreeOfFreedom(0).getEqId();

          
          this->homogenizationMatrix(posA, 7) = dx(2);  // x
          this->homogenizationMatrix(posA, 8) = -dx(1); // x
          this->homogenizationMatrix(posA, 9) = dx(0);  // x
          
          if (abs(dx(0)) <= prec(1e-10) && abs(dx(2)) > prec(1e-10)) {
            this->homogenizationMatrix(posA, 10) = coor(0) * dx(2); // x
          } else if (abs(dx(2)) < prec(1e-10)) {
            this->homogenizationMatrix(posA, 10) = dx(0) * coor(2); // x
          } else {
            this->homogenizationMatrix(posA, 10) = prec(0); // x
          }
          if (abs(dx(0)) > prec(1e-10)) {
            this->homogenizationMatrix(posA, 11) = -dx(0) * coor(1); // x
          } else {
            this->homogenizationMatrix(posA, 11) = prec(0); // x
          }

        }
      }
    }

    // Thermal periodic zEdges
    for (auto &sEdges : {zEdges2, zEdges3, zEdges4}) {
      for (auto i = 0; i < zEdges1.size(); ++i) {
        auto &mEdge = pointers.getGeometryData()->getEdgeData(zEdges1[i]);
        auto &sEdge1 = pointers.getGeometryData()->getEdgeData(sEdges[i]);

        for (auto i = 0; i < 2; ++i) {
          auto &mV = *mEdge.getVertex(i);
          auto &sV1 = *sEdge1.getVertex(i);
          Types::Vector3<prec> dx = sV1.getCoordinates() - mV.getCoordinates();
          Types::Vector3<prec> coor = sV1.getCoordinates();
          auto sNodes = sV1.getNodesOfSet(meshIdTemp);
          indexType posA = sNodes[0]->getDegreeOfFreedom(0).getEqId();

          
          this->homogenizationMatrix(posA, 7) = dx(2);  // x
          this->homogenizationMatrix(posA, 8) = -dx(1); // x
          this->homogenizationMatrix(posA, 9) = dx(0);  // x

          if (abs(dx(0)) <= prec(1e-10) && abs(dx(1)) > prec(1e-10)) {
            this->homogenizationMatrix(posA, 11) = -coor(0) * dx(1); // x
          } else if (abs(dx(1)) < prec(1e-10)) {
            this->homogenizationMatrix(posA, 11) = -dx(0) * coor(1); // x
          } else {
            this->homogenizationMatrix(posA, 11) = prec(0); // x
          }
          if (abs(dx(0)) > prec(1e-10)) {
            this->homogenizationMatrix(posA, 10) = dx(0) * coor(2); // x
          } else {
            this->homogenizationMatrix(posA, 10) = prec(0); // x
          }
          
        }
      }
    }

    // Thermal vertices
    for (auto i : {v1, v2, v3, v4, v5, v6, v7, v8}) {
      auto &V = pointers.getGeometryData()->getVertexData(i);
      auto nodes = V.getNodesOfSet(meshIdTemp);
      indexType posA = nodes[0]->getDegreeOfFreedom(0).getEqId();
      Types::Vector3<prec> coor = V.getCoordinates();
      this->homogenizationMatrix(posA, 6) = prec(1.0);           // x
      this->homogenizationMatrix(posA, 7) = coor(2);             // z
      this->homogenizationMatrix(posA, 8) = -coor(1);            // -y
      this->homogenizationMatrix(posA, 9) = coor(0);             // x
      this->homogenizationMatrix(posA, 10) = coor(0) * coor(2);  // x*z
      this->homogenizationMatrix(posA, 11) = -coor(0) * coor(1); // -x*y
      this->homogenizationMatrixDisp(posA, 6) = prec(1.0);           // x
      this->homogenizationMatrixDisp(posA, 7) = coor(2);             // z
      this->homogenizationMatrixDisp(posA, 8) = -coor(1);            // -y
      this->homogenizationMatrixDisp(posA, 9) = coor(0);             // x
      this->homogenizationMatrixDisp(posA, 10) = coor(0) * coor(2);  // x*z
      this->homogenizationMatrixDisp(posA, 11) = -coor(0) * coor(1); // -x*y
    }

    // Thermal edges x
  }

  pointers.getSPDLogger().trace("Computed A-Matrix:\n{}",
                                this->homogenizationMatrix);
}

auto Homogenization3DThermoMechBeam::getNumberOfStrains() -> indexType {
  return indexType(12);
}

auto Homogenization3DThermoMechBeam::getAMatrix() -> Types::MatrixXX<prec> & {
  return this->homogenizationMatrix;
}

auto Homogenization3DThermoMechBeam::getDv() -> prec {
  return (xmax(0) - xmin(0));
}

auto Homogenization3DThermoMechBeam::getDisplacementIncrement(
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
  } else if (bctype == 2) {
    disp = this->homogenizationMatrixDisp * strainIncrement;
  }

  return disp;
}

void Homogenization3DThermoMechBeam::setPeriodicDisplacements(
    PointerCollection &pointers, Types::VectorX<prec> &strains,
    Types::VectorX<prec> &strainIncrement) {

  auto setDb = [](PointerCollection &pointers, DegreeOfFreedom &Dof, prec B) {
    if (Dof.getStatus() == dofStatus::constraint) {
      auto bConst =
          pointers.getSolutionState()->getConstraint(Dof.getConstraintId());
      auto genLink = reinterpret_cast<std::shared_ptr<GeneralLink> &>(bConst);
      genLink->setB(B);
    }
  };
  if (bctype == 1||bctype==2) {
    Types::Vector3<prec> u;
    Types::Matrix3X<prec> A;
    A.resize(3, strains.rows());
    A.setZero();
    for (indexType i = 0; i < static_cast<indexType>(leftFacesMaster.size()); ++i) {
      auto mFace = pointers.getGeometryData()->getFaceData(leftFacesMaster[i]);
      auto sFace = pointers.getGeometryData()->getFaceData(rightFacesSlave[i]);

      indexType numV = mFace->getNumberOfVerts();
      for (indexType nn = 0; nn < numV; ++nn) {
        auto mVert = mFace->getVertex(nn);
        auto sVert = sFace->getVertex(nn);
        Types::Vector3<prec> sCoor = sVert->getCoordinates();
        Types::Vector3<prec> dx = sCoor - mVert->getCoordinates();
        auto sNodes = sVert->getNodesOfSet(meshIdDisp);
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
  // Add Temperature perdiodic bc
  if (bctype == 2) {
    // left right faces
    Types::VectorXT<prec> A;
    A.resize(12);
    A.setZero();
    for (indexType i = 0; i < static_cast<indexType>(leftFacesMaster.size()); ++i) {
      auto mFace = pointers.getGeometryData()->getFaceData(leftFacesMaster[i]);
      auto sFace = pointers.getGeometryData()->getFaceData(rightFacesSlave[i]);

      indexType numV = mFace->getNumberOfVerts();
      for (indexType nn = 0; nn < numV; ++nn) {
        auto mVert = mFace->getVertex(nn);
        auto sVert = sFace->getVertex(nn);
        Types::Vector3<prec> sCoor = sVert->getCoordinates();
        Types::Vector3<prec> dx = sCoor - mVert->getCoordinates();
        auto sNodes = sVert->getNodesOfSet(meshIdTemp);
        A(7) = dx(2);
        A(8) = -dx(1);
        A(9) = dx(0);
        A(10) = dx(0) * sCoor(2);
        A(11) = -dx(0) * sCoor(1);
        auto u = A * strains;
        auto &DofA = sNodes[0]->getDegreeOfFreedom(0);
        setDb(pointers, DofA, u(0));
      }
    }
    // top bottom faces
    A.setZero();
    for (indexType i = 0; i < static_cast<indexType>(bottomFacesMaster.size()); ++i) {
      auto mFace = pointers.getGeometryData()->getFaceData(bottomFacesMaster[i]);
      auto sFace = pointers.getGeometryData()->getFaceData(topFacesSlave[i]);

      indexType numV = mFace->getNumberOfVerts();
      for (indexType nn = 0; nn < numV; ++nn) {
        auto mVert = mFace->getVertex(nn);
        auto sVert = sFace->getVertex(nn);
        Types::Vector3<prec> sCoor = sVert->getCoordinates();
        Types::Vector3<prec> dx = sCoor - mVert->getCoordinates();
        auto sNodes = sVert->getNodesOfSet(meshIdTemp);
        A(7) = dx(2);
        A(10) = sCoor(0) * dx(2);
        auto u = A * strains;
        auto &DofA = sNodes[0]->getDegreeOfFreedom(0);
        setDb(pointers, DofA, u(0));
      }
    }
    // front back faces
    A.setZero();
    for (indexType i = 0; i < static_cast<indexType>(backFacesMaster.size()); ++i) {
      auto mFace = pointers.getGeometryData()->getFaceData(backFacesMaster[i]);
      auto sFace = pointers.getGeometryData()->getFaceData(frontFacesSlave[i]);

      indexType numV = mFace->getNumberOfVerts();
      for (indexType nn = 0; nn < numV; ++nn) {
        auto mVert = mFace->getVertex(nn);
        auto sVert = sFace->getVertex(nn);
        Types::Vector3<prec> sCoor = sVert->getCoordinates();
        Types::Vector3<prec> dx = sCoor - mVert->getCoordinates();
        auto sNodes = sVert->getNodesOfSet(meshIdTemp);
        A(8) = -dx(1);
        A(11) = -sCoor(0) * dx(1);
        auto u = A * strains;
        auto &DofA = sNodes[0]->getDegreeOfFreedom(0);
        setDb(pointers, DofA, u(0));
      }
    }


    // xEdges
    A.setZero();
    for (auto &sEdges : {xEdges2, xEdges3, xEdges4}) {
      for (auto i = 0; i < xEdges1.size(); ++i) {
        auto &mEdge = pointers.getGeometryData()->getEdgeData(xEdges1[i]);
        auto &sEdge1 = pointers.getGeometryData()->getEdgeData(sEdges[i]);

        for (auto i = 0; i < 2; ++i) {
          auto &mV = *mEdge.getVertex(i);
          auto &sV1 = *sEdge1.getVertex(i);
          auto dX = sV1.getCoordinates() - mV.getCoordinates();
          auto sCoor = sV1.getCoordinates();
          A(7) = dX(2);
          A(8) = -dX(1);
          A(9) = dX(0);
          A(10) = sCoor(0) * dX(2);
          A(11) = -sCoor(0) * dX(1);

          auto u = A * strains;
          auto Nodes = sV1.getNodesOfSet(meshIdTemp);
          auto &Dof = Nodes[0]->getDegreeOfFreedom(0);
          if (Dof.getStatus() == dofStatus::constraint) {
            setDb(pointers, Dof, u(0));
          }
        }
      }
    }
    // yEdges
    A.setZero();
    for (auto &sEdges : {yEdges2, yEdges3, yEdges4}) {
      for (auto i = 0; i < yEdges1.size(); ++i) {
        auto &mEdge = pointers.getGeometryData()->getEdgeData(yEdges1[i]);
        auto &sEdge1 = pointers.getGeometryData()->getEdgeData(sEdges[i]);

        for (auto i = 0; i < 2; ++i) {
          auto &mV = *mEdge.getVertex(i);
          auto &sV1 = *sEdge1.getVertex(i);
          auto dX = sV1.getCoordinates() - mV.getCoordinates();
          auto sCoor = sV1.getCoordinates();
          A(7) = dX(2);
          A(8) = -dX(1);
          A(9) = dX(0);
          if (abs(dX(0)) <= prec(1e-10) && abs(dX(2)) > prec(1e-10)) {
            A(10) = sCoor(0) * dX(2);
          } else if (abs(dX(2)) < prec(1e-10)) {
            A(10) = dX(0) * sCoor(2);
          } else
          {
            A(10) = prec(0);
          }
          if (abs(dX(0))>prec(1e-10))
          {
            A(11) = -dX(0) * sCoor(1);
          } else
          {
            A(11) = prec(0);
          }
          auto u = A * strains;
          auto Nodes = sV1.getNodesOfSet(meshIdTemp);
          auto &Dof = Nodes[0]->getDegreeOfFreedom(0);
          if (Dof.getStatus() == dofStatus::constraint) {
            setDb(pointers, Dof, u(0));
          }
        }
      }
    }
    // z edges
    A.setZero();
    for (auto &sEdges : {zEdges2, zEdges3, zEdges4}) {
      for (auto i = 0; i < zEdges1.size(); ++i) {
        auto &mEdge = pointers.getGeometryData()->getEdgeData(zEdges1[i]);
        auto &sEdge1 = pointers.getGeometryData()->getEdgeData(sEdges[i]);

        for (auto i = 0; i < 2; ++i) {
          auto &mV = *mEdge.getVertex(i);
          auto &sV1 = *sEdge1.getVertex(i);
          auto dX = sV1.getCoordinates() - mV.getCoordinates();
          auto sCoor = sV1.getCoordinates();
          A(7) = dX(2);
          A(8) = -dX(1);
          A(9) = dX(0);
          if (abs(dX(0)) <= prec(1e-10) && abs(dX(1)) > prec(1e-10)) {
            A(11) = -sCoor(0) * dX(1);
          } else if (abs(dX(1)) < prec(1e-10)) {
            A(11) = -dX(0) * sCoor(1);
          } else {
            A(11) = prec(0);
          }
          if (abs(dX(0)) > prec(1e-10)) 
          {
            A(10) = dX(0) * sCoor(2);
          } else
          {
            A(10) = prec(0);
          }
          auto u = A * strains;
          auto Nodes = sV1.getNodesOfSet(meshIdTemp);
          auto &Dof = Nodes[0]->getDegreeOfFreedom(0);
          if (Dof.getStatus() == dofStatus::constraint) {
            setDb(pointers, Dof, u(0));
          }
        }
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
  writeEigenMatrix(out, xmin);
  writeEigenMatrix(out, xmax);
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
  readEigenMatrix(in, xmin);
  readEigenMatrix(in, xmax);
}

} // namespace HierAMuS