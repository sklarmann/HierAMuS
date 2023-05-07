// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause




#include "plot/vtkplotClassBase.h"
#include <plot/vtkplotClass.h>


#include <types/MatrixTypes.h>

#include <elementFormulations/GenericElementFormulation.h>
#include <equations/NodeSet.h>
#include <finiteElements/LinearPrism.h>
#include <geometry/Base.h>
#include <geometry/Edges.h>
#include <geometry/GeometryData.h>
#include <geometry/Vertex.h>
#include <geometry/Faces.h>
#include <geometry/Volumes.h>
#include <solver/GenericSolutionState.h>

#include <pointercollection/pointercollection.h>

#include "shapefunctions/IntegrationsPoints/dataClasses/GaussPoints.h"
#include <shapefunctions/KernelShapes.h>

#include "vtkCellType.h"



namespace HierAMuS::FiniteElement {
  

LinearPrism::~LinearPrism() = default;

auto LinearPrism::getVertexIds(PointerCollection& pointers) -> std::vector<indexType> {

  std::vector<indexType> u;
  Geometry::Base *geoElem;
  geoElem = pointers.getGeometryData()->getVolume(this->vol);
  geoElem->getVerts(u);
  return u;
}

void LinearPrism::setVolume(indexType volIn) { this->vol = volIn; }

auto LinearPrism::getVertex(ptrCol &pointers, indexType num)
    -> Geometry::Vertex & {

  std::vector<indexType> vertIds;
  vertIds = this->getVertexIds(pointers);

  return pointers.getGeometryData()->getVertex(vertIds[num]);
}

auto LinearPrism::getEdge(ptrCol &pointers, indexType num)
    -> Geometry::Edges & {
  Geometry::Base *temp;
  temp = pointers.getGeometryData()->getVolume(this->vol);
  std::vector<indexType> EdgeNums;
  temp->getEdges(EdgeNums);
  return pointers.getGeometryData()->getEdge(EdgeNums[num]);
}


void LinearPrism::getH1Shapes(ptrCol &pointers, indexType order,
                              Types::Matrix33<prec> &jacobi,
                              Types::VectorX<prec> &shape,
                              Types::Matrix3X<prec> &dshape, prec xi,
                              prec eta, prec zeta) {
  Geometry::Volumes *tempVol;
  tempVol = pointers.getGeometryData()->getVolume(this->vol);
  tempVol->getH1Shapes(pointers, order, shape, dshape, xi, eta, zeta);
  dshape = jacobi.inverse().transpose() * dshape;
}

void LinearPrism::getJacobian(ptrCol &pointers, Types::Matrix33<prec> &jacobi,
                              prec xsi, prec eta,
                              prec zeta) {
  Geometry::Volumes *tempVol;
  tempVol = pointers.getGeometryData()->getVolume(this->vol);
  tempVol->getJacobian(pointers, jacobi, xsi, eta, zeta);
}

void LinearPrism::setH1Shapes(ptrCol &pointers, indexType meshid,
                              indexType order) {
  Geometry::Volumes *tempVol;
  tempVol = pointers.getGeometryData()->getVolume(this->vol);
  tempVol->setH1Shapes(pointers, meshid, order, NodeTypes::displacement);
}

void LinearPrism::getH1Dofs(ptrCol &pointers,
                            std::vector<DegreeOfFreedom *> &Dofs,
                            indexType meshID, indexType order) {
  Geometry::Volumes *temp = pointers.getGeometryData()->getVolume(this->vol);
  temp->getH1Dofs(pointers, Dofs, meshID, order);
}

void LinearPrism::getGaussPoints(indexType number,
                                        std::vector<prec> &weight,
                                        std::vector<prec> &xsi,
                                        std::vector<prec> &eta,
                                        std::vector<prec> &zeta) {
  brickGP(xsi, eta, zeta, weight, number);
}



void LinearPrism::toParaviewAdapter(PointerCollection &ptrCol,
                                    vtkPlotInterface &catalyst,
                                    const ParaviewSwitch &ToDo) {
#ifdef USE_VTK

  switch (ToDo) {
  case ParaviewSwitch::Mesh: {
    int matNum = static_cast<int>(this->getMaterial()->getNumber());
    Types::Vector3<prec> coors(3);
    std::vector<double> coorsA(3);
    Geometry::Base *Vert;
    std::vector<indexType> cell(6);
    for (indexType i = 0; i < 6; ++i) {
      Vert = &this->getVertex(ptrCol, i);
      coors = Vert->getCoordinates();
      for (auto j = 0; j < 3; ++j) {
        coorsA[j] = static_cast<double>(coors(j));
      }
      catalyst.addPoint(0, matNum, Vert->getId(), coorsA[0], coorsA[1],
                        coorsA[2]);
      cell[i] = Vert->getId();
    }
    int celltype = VTK_WEDGE;
    indexType dummy = 1;
    indexType nn = 6;
    catalyst.addCell(0, matNum, this->id, dummy, cell, nn, celltype);
    break;
  }
  case ParaviewSwitch::Solution: {
    indexType matNum = (this->getMaterial()->getNumber());
    Geometry::Base *Vert;
    NodeSet *tempSet;
    std::vector<DegreeOfFreedom *> Dofs;
    Types::VectorX<prec> solution;
    std::vector<prec> sol(3);
    for (indexType i = 0; i < 6; ++i) {
      Vert = &this->getVertex(ptrCol, i);
      std::vector<NodeSet *> sets;
      Vert->getSets(ptrCol, sets);
      indexType vertId = (Vert->getId());
      for (auto j = sets.begin(); j != sets.end(); ++j) {
        tempSet = *j;
        indexType id = tempSet->getMeshId();
        indexType nnodes = tempSet->getNumberOfNodes();
        for (auto k = 0; k < nnodes; ++k) {
          std::stringstream ArrName;
          ArrName << "SolutionM" << id << "N" << k;
          Dofs = tempSet->getDegreesOfFreedom(ptrCol);
          this->getSolution(ptrCol, Dofs, solution);
          for (auto i = 0; i < 3; ++i) {
            sol[i] = static_cast<double>(solution(i));
          }
          indexType comp = 3;
          indexType main = 0;
          catalyst.setPointData(main, matNum, vertId, sol, comp, ArrName.str());
        }
      }
    }
    break;
  }
  case ParaviewSwitch::Weights: {
  } break;
  case ParaviewSwitch::ProjectedValues: {
  } break;
  case ParaviewSwitch::Eigenvector:
    break;
  }

  // select(ToDo){
  //   case ParaviewSwitch::Mesh:
  //   int matNum = static_cast<int>(this->getMaterial()->getNumber());
  //   std::vector<prec> coors;
  //   GenericGeometryElement *Vert;
  //   std::vector<int> cell(8);
  //   for(indexType i=0;i<8;++i){
  //     Vert = this->getVertex(ptrCol,i);
  //     coors = Vert->getCoordinates();
  //     int id = static_cast<int>(Vert->getId());
  //     catalyst.addPoint(0,matNum,id,coors[0],coors[1],coors[2]);
  //     cell[i] = static_cast<int>(id);
  //   }
  //   int celltype = VTK_HEXAHEDRON;
  //   int dummy = 1;
  //   int nn = 8;
  //   int idd = this->id;
  //   catalyst.addCell(0,matNum,idd,dummy,&cell[0],&nn,&celltype);
  //   break;

  //   case ParaviewSwitch::Solution:
  //   break;
  // }

#endif
}

void LinearPrism::projectOnVertsParaview(PointerCollection &ptrCol,
                                         vtkPlotInterface &catalyst,
                                         Types::VectorX<prec> &values,
                                         prec &xsi, prec &eta,

                                         prec &zeta, prec &weight,
                                         std::string name) {
  Types::Matrix3X<prec> shapeDeriv;
  Types::VectorX<prec> shape;
  Types::Matrix33<prec> jacobi;
  this->getJacobian(ptrCol, jacobi, xsi, eta, zeta);
  this->getH1Shapes(ptrCol, 1, jacobi, shape, shapeDeriv, xsi, eta, zeta);
  prec detj = jacobi.determinant();
  prec dvp;
  std::vector<indexType> verts = this->getVertexIds(ptrCol);
  std::vector<prec> toAdd(values.size());
  indexType numComp = values.size();
  for (auto i = 0; i < verts.size(); i++) {
    dvp = detj * weight * shape(i);
    for (auto i = 0; i < values.size(); ++i) {
      toAdd[i] = values(i) * dvp;
    }
    catalyst.SumPointDataWeighted(0, this->getMaterial()->getNumber(), toAdd,
                                  verts[i], numComp, name);
  }
}

void LinearPrism::setH1BeamShapes(LinearPrism::ptrCol &pointers,
                                  indexType meshid, indexType order) {

  auto &edge = this->getBeamEdge(pointers);

  edge.setH1Shapes(pointers, meshid, order, NodeTypes::displacement);
}

void LinearPrism::getH1BeamDofs(LinearPrism::ptrCol &pointers,
                                std::vector<DegreeOfFreedom *> &Dofs,
                                indexType meshID, indexType order) {
  auto &edge = this->getBeamEdge(pointers);

  edge.getH1Dofs(pointers, Dofs, meshID, order);
}

void LinearPrism::getH1BeamShapes(LinearPrism::ptrCol &pointers,
                                  indexType order, prec jacobian,
                                  Types::VectorX<prec> &shape,
                                  Types::VectorX<prec> &shapeDerivatives,
                                  prec xi) {

  auto &edge = this->getBeamEdge(pointers);

  edge.getH1Shapes(pointers, order, shape, shapeDerivatives, xi);
  shapeDerivatives /= jacobian;
}

auto LinearPrism::getBeamEdge(LinearPrism::ptrCol &pointers)
    -> Geometry::Edges & {
  Geometry::Volumes *geoElem = pointers.getGeometryData()->getVolume(this->vol);

  std::vector<indexType> edges;
  geoElem->getEdges(edges);

  auto &edge = pointers.getGeometryData()->getEdge(edges[3]);
  return edge;
}

auto LinearPrism::getEndTriad(LinearPrism::ptrCol &pointers)
    -> Types::Matrix33<prec> {

  Geometry::Volumes *prism = pointers.getGeometryData()->getVolume(this->vol);

  std::vector<indexType> faceNums;
  prism->getFaces(faceNums);

  Geometry::Faces *triad;
  triad = pointers.getGeometryData()->getFace(faceNums.back());

  faceNums.clear();

  triad->getVerts(faceNums);



  auto &v1 = pointers.getGeometryData()->getVertex(faceNums[0]);
  auto &v2 = pointers.getGeometryData()->getVertex(faceNums[1]);
  auto &v3 = pointers.getGeometryData()->getVertex(faceNums[2]);

  Types::Vector3<prec> A2 = v2.getCoordinates() - v1.getCoordinates();
  Types::Vector3<prec> A3 = v3.getCoordinates() - v1.getCoordinates();
  Types::Vector3<prec> A1 = A2.cross(A3);
  A1 = A1.normalized();
  A2 = A2.normalized();
  A3 = A3.normalized();
  Types::Matrix33<prec> retMatrix;

  retMatrix.block(0, 0, 3, 1) = A1;
  retMatrix.block(0, 1, 3, 1) = A2;
  retMatrix.block(0, 2, 3, 1) = A3;

  return retMatrix;
}

auto LinearPrism::getStartTriad(LinearPrism::ptrCol &pointers)
    -> Types::Matrix33<prec> {
  Geometry::Volumes *prism = pointers.getGeometryData()->getVolume(this->vol);

  std::vector<indexType> faceNums;
  prism->getFaces(faceNums);

  Geometry::Faces *triad;
  triad = pointers.getGeometryData()->getFace(faceNums.front());

  faceNums.clear();

  triad->getVerts(faceNums);


  auto &v1 = pointers.getGeometryData()->getVertex(faceNums[0]);
  auto &v2 = pointers.getGeometryData()->getVertex(faceNums[1]);
  auto &v3 = pointers.getGeometryData()->getVertex(faceNums[2]);

  Types::Vector3<prec> A2 = v2.getCoordinates() - v1.getCoordinates();
  Types::Vector3<prec> A3 = v3.getCoordinates() - v1.getCoordinates();
  Types::Vector3<prec> A1 = A2.cross(A3);
  A1 = A1.normalized();
  A2 = A2.normalized();
  A3 = A3.normalized();
  Types::Matrix33<prec> retMatrix;

  retMatrix.block(0, 0, 3, 1) = A1;
  retMatrix.block(0, 1, 3, 1) = A2;
  retMatrix.block(0, 2, 3, 1) = A3;

  return retMatrix;
}

auto LinearPrism::getBeamJacobian(LinearPrism::ptrCol &pointers, prec xi)
    -> prec {

  prec jac;

  auto &edge = this->getBeamEdge(pointers);

  jac = edge.getJacobian(pointers, xi);

  return jac;
};

} // namespace HierAMuS::FiniteElement

