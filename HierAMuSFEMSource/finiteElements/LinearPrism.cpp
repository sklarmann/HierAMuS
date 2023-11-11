// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "geometry/GeometryData.h"


#include "plot/vtkplotClassBase.h"
#include <plot/vtkplotClass.h>


#include <types/MatrixTypes.h>

#include <elementFormulations/GenericElementFormulation.h>
#include <finiteElements/LinearPrism.h>
#include <geometry/GeometryBaseData.h>
#include <geometry/Edges/EdgesData.h>
#include <geometry/VertexData.h>
#include <geometry/Faces/FacesData.h>
#include <geometry/Volumes/VolumesData.h>
#include <solver/GenericSolutionState.h>

#include "geometry/Volumes/VolumesRuntime.h"
#include "geometry/Volumes/VolumesH1Interface.h"

#include <pointercollection/pointercollection.h>

#include "shapefunctions/IntegrationsPoints/dataClasses/GaussPoints.h"
#include <shapefunctions/KernelShapes.h>

#include "vtkCellType.h"



namespace HierAMuS::FiniteElement {
  

LinearPrism::~LinearPrism() = default;

void LinearPrism::set_pointers(PointerCollection &pointers) {
  m_volume_runtime =
      pointers.getGeometryData()->getVolumeRuntime(m_volume);
}

auto LinearPrism::getVertexIds(PointerCollection& pointers) -> std::vector<indexType> {

  auto geoElem = pointers.getGeometryData()->getVolumeData(this->m_volume);
  return geoElem->getVertexNumbers();
}

void LinearPrism::setVolume(indexType volIn) { this->m_volume = volIn; }

auto LinearPrism::getVertex(ptrCol &pointers, indexType num)
    -> Geometry::VertexData & {

  std::vector<indexType> vertIds;
  vertIds = this->getVertexIds(pointers);

  return pointers.getGeometryData()->getVertexData(vertIds[num]);
}

auto LinearPrism::getEdge(ptrCol &pointers, indexType num)
    -> Geometry::EdgesData & {
  auto temp = pointers.getGeometryData()->getVolumeData(this->m_volume);
  std::vector<indexType> EdgeNums;
  temp->getEdgeNumbers(EdgeNums);
  return pointers.getGeometryData()->getEdgeData(EdgeNums[num]);
}


void LinearPrism::getH1Shapes(ptrCol &pointers, indexType order,
                              Types::Matrix33<prec> &jacobi,
                              Types::VectorX<prec> &shape,
                              Types::Matrix3X<prec> &dshape, prec xi,
                              prec eta, prec zeta) {
  m_volume_runtime->getH1Volume()->getH1Shapes(order, shape, dshape,
                                               xi, eta, zeta);
  dshape = jacobi.inverse().transpose() * dshape;
}

void LinearPrism::getJacobian(ptrCol &pointers, Types::Matrix33<prec> &jacobi,
                              prec xsi, prec eta,
                              prec zeta) {
  IntegrationPoint a;
  a.xi = xsi;
  a.eta = eta;
  a.zeta = zeta;
  jacobi = m_volume_runtime->getJacobian(a);
}

void LinearPrism::setH1Shapes(ptrCol &pointers, indexType meshid,
                              indexType order) {
  Geometry::VolumesData *tempVol;
  tempVol = pointers.getGeometryData()->getVolumeData(this->m_volume);
  tempVol->setH1Shapes(meshid, order, NodeTypes::displacement);
}

void LinearPrism::getH1Dofs(ptrCol &pointers,
                            std::vector<DegreeOfFreedom *> &Dofs,
                            indexType meshID, indexType order) {
  m_volume_runtime->getH1Volume()->getH1Dofs(Dofs, meshID, order);
}

void LinearPrism::getGaussPoints(indexType number,
                                        std::vector<prec> &weight,
                                        std::vector<prec> &xsi,
                                        std::vector<prec> &eta,
                                        std::vector<prec> &zeta) {
  brickGP(xsi, eta, zeta, weight, number);
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

  edge.setH1Shapes(meshid, order, NodeTypes::displacement);
}

void LinearPrism::getH1BeamDofs(LinearPrism::ptrCol &pointers,
                                std::vector<DegreeOfFreedom *> &Dofs,
                                indexType meshID, indexType order) {
  auto &edge = this->getBeamEdge(pointers);

  edge.getH1Dofs(Dofs, meshID, order);
}

void LinearPrism::getH1BeamShapes(LinearPrism::ptrCol &pointers,
                                  indexType order, prec jacobian,
                                  Types::VectorX<prec> &shape,
                                  Types::VectorX<prec> &shapeDerivatives,
                                  prec xi) {

  auto &edge = this->getBeamEdge(pointers);

  edge.getH1Shapes(order, shape, shapeDerivatives, xi);
  shapeDerivatives /= jacobian;
}

auto LinearPrism::getBeamEdge(LinearPrism::ptrCol &pointers)
    -> Geometry::EdgesData & {
  Geometry::VolumesData *geoElem = pointers.getGeometryData()->getVolumeData(this->m_volume);

  std::vector<indexType> edges;
  geoElem->getEdgeNumbers(edges);

  auto &edge = pointers.getGeometryData()->getEdgeData(edges[3]);
  return edge;
}

auto LinearPrism::getEndTriad(LinearPrism::ptrCol &pointers)
    -> Types::Matrix33<prec> {

  Geometry::VolumesData *prism = pointers.getGeometryData()->getVolumeData(this->m_volume);

  std::vector<indexType> faceNums;
  prism->getFaceNumbers(faceNums);

  Geometry::FacesData *triad;
  triad = pointers.getGeometryData()->getFaceData(faceNums.back());

  faceNums.clear();

  faceNums = triad->getVertexNumbers();



  auto &v1 = pointers.getGeometryData()->getVertexData(faceNums[0]);
  auto &v2 = pointers.getGeometryData()->getVertexData(faceNums[1]);
  auto &v3 = pointers.getGeometryData()->getVertexData(faceNums[2]);

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
  Geometry::VolumesData *prism = pointers.getGeometryData()->getVolumeData(this->m_volume);

  std::vector<indexType> faceNums;
  prism->getFaceNumbers(faceNums);

  Geometry::FacesData *triad;
  triad = pointers.getGeometryData()->getFaceData(faceNums.front());

  faceNums.clear();

  faceNums = triad->getVertexNumbers();


  auto &v1 = pointers.getGeometryData()->getVertexData(faceNums[0]);
  auto &v2 = pointers.getGeometryData()->getVertexData(faceNums[1]);
  auto &v3 = pointers.getGeometryData()->getVertexData(faceNums[2]);

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

  jac = edge.getJacobian(xi);

  return jac;
};

auto LinearPrism::getType() -> Elementtypes {
  return Elementtypes::LinearPrism;
}

} // namespace HierAMuS::FiniteElement

