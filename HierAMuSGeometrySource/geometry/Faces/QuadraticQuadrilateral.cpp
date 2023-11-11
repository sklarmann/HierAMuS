// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "MatrixTypes.h"
#include "datatypes.h"
#include "geometry/GeometryBaseData.h"
#include "geometry/Faces/FacesData.h"
#include "geometry/GeometryTypes.h"
#include "geometry/Faces/LinearQuadrilateralData.h"
#include "plot/vtkplotClassBase.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "shapefunctions/LagrangeShape.h"
#include "shapefunctions/LegendreShapes.h"
#include "shapefunctions/LobattoShapes.h"
#include "exception"

#include "control/HandlingStructs.h"
#include "control/OutputHandler.h"
#include "geometry/Faces/QuadraticQuadrilateral.h"
#include "pointercollection/pointercollection.h"
#include "sstream"

#include "plot/vtkplotClass.h"

#include "geometry/GeometryData.h"
#include "geometry/Edges/EdgesData.h"
#include "geometry/GeometryData.h"
#include "geometry/VertexData.h"

#include "LoadList.h"
#include "solver/GenericSolutionState.h"

#include "stdexcept"
#include "vector"
#include "vtkCellType.h"

#include "HelperFunctions.h"

namespace HierAMuS::Geometry {

QuadraticQuadrilateral::QuadraticQuadrilateral()
    : FacesDataInterface(){};

QuadraticQuadrilateral::~QuadraticQuadrilateral() = default;

auto QuadraticQuadrilateral::getType() -> const GeometryTypes & {
  return this->type;
}


auto QuadraticQuadrilateral::getCoordinates(PointerCollection &pointers,
                                            prec xi, prec eta)
    -> Types::Vector3<prec> {
  Types::Vector3<prec> ret;
  Types::Vector3<prec> t1;
  ret.setZero();
  std::array<prec, 3> eta_s;
  std::array<prec, 3> eta_ds;
  std::array<prec, 3> xi_s;
  std::array<prec, 3> xi_ds;
  std::array<prec, 9> vals;

  LagrangeShape(xi_s[1], xi_ds[1], xi, 2, 0);
  LagrangeShape(xi_s[2], xi_ds[2], xi, 2, 2);
  LagrangeShape(xi_s[3], xi_ds[3], xi, 2, 1);
  LagrangeShape(eta_s[1], eta_ds[1], eta, 2, 0);
  LagrangeShape(eta_s[2], eta_ds[2], eta, 2, 2);
  LagrangeShape(eta_s[3], eta_ds[3], eta, 2, 1);

  vals[0] = eta_s[1] * xi_s[1];
  vals[1] = eta_s[3] * xi_s[1];
  vals[2] = eta_s[3] * xi_s[3];
  vals[3] = eta_s[1] * xi_s[3];
  vals[4] = eta_s[2] * xi_s[1];
  vals[5] = eta_s[3] * xi_s[2];
  vals[6] = eta_s[2] * xi_s[3];
  vals[7] = eta_s[1] * xi_s[2];
  vals[8] = eta_s[2] * xi_s[2];

  for (auto i = 0; i < 4; ++i) {
    t1 = pointers.getGeometryData()->getVertexData(this->m_verts[i]).getCoordinates();
    ret += t1 * vals[i];
  }
  for (auto i = 4; i < 9; ++i) {
    t1 = pointers.getGeometryData()
             ->getVertexData(this->m_verts[i])
             .getCoordinates();
    ret += t1 * vals[i];
  }
  return ret;
}

auto QuadraticQuadrilateral::getCoordinates(PointerCollection &pointers,
                                            IntegrationPoint &integrationPoint)
    -> Types::Vector3<prec> {
  Types::Vector3<prec> ret;
  Types::Vector3<prec> t1;
  ret.setZero();
  std::array<prec, 3> eta_s;
  std::array<prec, 3> eta_ds;
  std::array<prec, 3> xi_s;
  std::array<prec, 3> xi_ds;
  std::array<prec, 9> vals;

  LagrangeShape(xi_s[1], xi_ds[1], integrationPoint.xi, 2, 0);
  LagrangeShape(xi_s[2], xi_ds[2], integrationPoint.xi, 2, 1);
  LagrangeShape(xi_s[3], xi_ds[3], integrationPoint.xi, 2, 2);
  LagrangeShape(eta_s[1], eta_ds[1], integrationPoint.eta, 2, 0);
  LagrangeShape(eta_s[2], eta_ds[2], integrationPoint.eta, 2, 1);
  LagrangeShape(eta_s[3], eta_ds[3], integrationPoint.eta, 2, 2);

  vals[0] = eta_s[1] * xi_s[1];
  vals[1] = eta_s[3] * xi_s[1];
  vals[2] = eta_s[3] * xi_s[3];
  vals[3] = eta_s[1] * xi_s[3];
  vals[4] = eta_s[2] * xi_s[1];
  vals[5] = eta_s[3] * xi_s[2];
  vals[6] = eta_s[2] * xi_s[3];
  vals[7] = eta_s[1] * xi_s[2];
  vals[8] = eta_s[2] * xi_s[2];

  for (auto i = 0; i < 8; ++i) {
    t1 = pointers.getGeometryData()->getVertexData(this->m_verts[i]).getCoordinates();
    ret += t1 * vals[i];
  }
  return ret;
}

auto QuadraticQuadrilateral::getJacobian(PointerCollection &pointers,
                                         IntegrationPoint &IntegrationPt)
    -> Types::Matrix22<prec> {
  Types::Matrix22<prec> jacobi;
  jacobi.setZero();

  std::array<prec, 3> eta_s;
  std::array<prec, 3> eta_ds;
  std::array<prec, 3> xi_s;
  std::array<prec, 3> xi_ds;
  std::array<prec, 9> N_de;
  std::array<prec, 9> N_dxi;

  LagrangeShape(xi_s[0], xi_ds[0], IntegrationPt.xi, 2, 0);
  LagrangeShape(xi_s[1], xi_ds[1], IntegrationPt.xi, 2, 2);
  LagrangeShape(xi_s[2], xi_ds[2], IntegrationPt.xi, 2, 1);
  LagrangeShape(eta_s[0], eta_ds[0], IntegrationPt.eta, 2, 0);
  LagrangeShape(eta_s[1], eta_ds[1], IntegrationPt.eta, 2, 2);
  LagrangeShape(eta_s[2], eta_ds[2], IntegrationPt.eta, 2, 1);

  N_de[0] = eta_ds[0] * xi_s[0];
  N_de[1] = eta_ds[0] * xi_s[2];
  N_de[2] = eta_ds[2] * xi_s[2];
  N_de[3] = eta_ds[2] * xi_s[0];
  N_de[4] = eta_ds[0] * xi_s[1];
  N_de[5] = eta_ds[1] * xi_s[2];
  N_de[6] = eta_ds[2] * xi_s[1];
  N_de[7] = eta_ds[1] * xi_s[0];
  N_de[8] = eta_ds[1] * xi_s[1];

  N_dxi[0] = eta_s[0] * xi_ds[0];
  N_dxi[1] = eta_s[0] * xi_ds[2];
  N_dxi[2] = eta_s[2] * xi_ds[2];
  N_dxi[3] = eta_s[2] * xi_ds[0];
  N_dxi[4] = eta_s[0] * xi_ds[1];
  N_dxi[5] = eta_s[1] * xi_ds[2];
  N_dxi[6] = eta_s[2] * xi_ds[1];
  N_dxi[7] = eta_s[1] * xi_ds[0];
  N_dxi[8] = eta_s[1] * xi_ds[1];

  for (auto i = 0; i < 4; ++i) {
    auto coord =
        pointers.getGeometryData()->getVertexData(this->m_verts[i]).getCoordinates();
    jacobi(0, 0) += N_dxi[i] * coord(0);
    jacobi(0, 1) += N_de[i] * coord(0);
    jacobi(1, 0) += N_dxi[i] * coord(1);
    jacobi(1, 1) += N_de[i] * coord(1);
  }
  for (auto i = 4; i < 9; ++i) {
    auto coord = pointers.getGeometryData()
                     ->getVertexData(this->m_verts[i])
                     .getCoordinates();
    jacobi(0, 0) += N_dxi[i] * coord(0);
    jacobi(0, 1) += N_de[i] * coord(0);
    jacobi(1, 0) += N_dxi[i] * coord(1);
    jacobi(1, 1) += N_de[i] * coord(1);
  }
  // jacobi = jacobi.transpose().eval();
  return jacobi;
}



void QuadraticQuadrilateral::geometryToParaview(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh) {
  indexType numPoints = 9;
  std::vector<indexType> points(numPoints);
  points.clear();
  // for (auto i : this->verts) {
  //   auto &vert = pointers.getGeometryData()->getVertex(i);
  //   vert.geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
  //   points.push_back(i);
  // }
  // for (auto i : this->verts_quad) {
  //   auto &vert = pointers.getGeometryData()->getVertex(i);
  //   vert.geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
  //   points.push_back(i);
  // }

  for (auto i = 0; i < 9; ++i) {
    auto &vert = *this->getVertex(pointers, i);
    vert.geometryToParaview(pointers, paraviewAdapter, mainMesh, subMesh);
    points.push_back(vert.getId());
  }

  paraviewAdapter.addCell(mainMesh, subMesh, this->id, 1, points, numPoints,
                          VTK_BIQUADRATIC_QUAD);
}

void QuadraticQuadrilateral::computeWeightsParaview(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh) {

  auto GP = this->getIntegrationPoints(-1);
  GP.setOrder(3);
  //auto N1 = [](prec xi) { return (xi - prec(1)) * xi / prec(2); };
  //auto N3 = [](prec xi) { return (prec(1) - xi) * (prec(1) + xi); };
  //auto N2 = [](prec xi) { return (xi + prec(1)) * xi / prec(2); };

  for (auto i : GP) {
    auto jaco = this->getJacobian(pointers, i);
    auto geoShapes = this->getGeometryShapes(pointers, i);
    prec dA = jaco.determinant() * i.weight;
    for (auto nn = 0; nn < 9; ++nn) {
      auto &V = *this->getVertex(pointers, nn);
      std::vector<prec> val;
      val.push_back(geoShapes.shapes(nn) * dA);
      paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, val, V.getId(), 1,
                                           paraviewNames::weightName());
    }

    /*auto jaco = this->getJacobian(pointers, i);
    prec dA = jaco.determinant() * i.weight;
    auto shpXi = shp(i.xi);
    auto shpEta = shp(i.eta);

    {
      std::array<indexType, 4> idxXi = {0, 1, 1, 0};
      std::array<indexType, 4> idxEta = {0, 0, 1, 1};
      for (auto j = 0; j < 4; ++j) {
        auto &vert = pointers.getGeometryData()->getVertex(this->verts[j]);
        std::vector<prec> val;
        prec shpval = shpXi[idxXi[j]] * shpEta[idxEta[j]] * dA;
        val.push_back(shpval);
        paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, val,
                                             vert.getId(), 1,
                                             paraviewNames::weightName());
      }
    }
    {
      std::array<indexType, 5> idxXi = {2, 1, 2, 0, 2};
      std::array<indexType, 5> idxEta = {0, 2, 1, 2, 2};
      for (auto j = 0; j < 5; ++j) {
        auto &vert = pointers.getGeometryData()->getVertex(this->verts_quad[j]);
        std::vector<prec> val;
        prec shpval = shpXi[idxXi[j]] * shpEta[idxEta[j]] * dA;
        val.push_back(shpval);
        paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, val,
                                             vert.getId(), 1,
                                             paraviewNames::weightName());
      }
    }*/
  }
}

void QuadraticQuadrilateral::H1SolutionToParaview(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType meshId, indexType order,
    std::string &name) {

  auto Dofs = this->getH1Dofs(pointers, meshId, order);
  auto solution = pointers.getSolutionState()->getSolution(Dofs);

  for (auto i = 0; i < 4; ++i) {
    auto &V = pointers.getGeometryData()->getVertexData(this->m_verts[i]);
    std::vector<prec> sol(3);
    for (auto j = 0; j < 3; ++j) {
      sol[j] = solution(3 * i + j);
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, V.getId(), sol, 3, name);
  }
  std::vector<prec> xi = {0.0, 1.0, 0.0, -1.0, 0.0};
  std::vector<prec> eta = {-1.0, 0.0, 1.0, 0.0, 0.0};
  for (auto i = 4; i < 9; ++i) {
    IntegrationPoint ip;
    ip.xi = xi[i - 4];
    ip.eta = eta[i - 4];
    auto shapes = this->getH1Shapes(pointers, order, ip);
    auto &V = pointers.getGeometryData()->getVertexData(this->m_verts[i]);
    std::vector<prec> sol = {0.0, 0.0, 0.0};
    for (auto j = 0; j < shapes.shapes.rows(); ++j) {
      for (auto k = 0; k < 3; ++k) {
        sol[k] += shapes.shapes(j) * solution(3 * j + k);
      }
    }
    paraviewAdapter.setPointData(mainMesh, subMesh, V.getId(), sol, 3, name);
  }
}

void QuadraticQuadrilateral::H1DataToParaview(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, Types::VectorX<prec> &Data,
    indexType numberComponents, indexType order, std::string &name) {

  std::array<prec, 9> xi = {-1.0, 1.0, 1.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0};
  std::array<prec, 9> eta = {-1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 1.0, 0.0, 0.0};
  std::vector<prec> interData(numberComponents);

  for (auto i = 0; i < 9; ++i) {
    IntegrationPoint ip;
    ip.xi = xi[i];
    ip.eta = eta[i];
    auto shapes = this->getH1Shapes(pointers, order, ip);
    for (auto j = 0; j < numberComponents; ++j) {
      interData[j] = prec(0);
    }
    for (auto j = 0; j < shapes.shapes.rows(); ++j) {
      for (auto k = 0; k < numberComponents; ++k) {
        interData[k] += shapes.shapes(j) * Data(numberComponents * j + k);
      }
    }
    auto &Vert = *this->getVertex(pointers, i);

    paraviewAdapter.setPointData(mainMesh, subMesh, Vert.getId(), interData,
                                 numberComponents, name);
  }
}

void QuadraticQuadrilateral::projectDataToParaviewVertices(
    PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
    indexType mainMesh, indexType subMesh, indexType order,
    IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
    indexType numberComponents, std::string name) {

  auto geoShapes = this->getGeometryShapes(pointers, IntegrationPt);
  auto jaco = this->getJacobian(pointers, IntegrationPt);
  auto dA = jaco.determinant() * IntegrationPt.weight;
  std::vector<prec> vals(numberComponents);
  for (auto i = 0; i < 9; ++i) {
    auto &V = *this->getVertex(pointers, i);
    for (auto j = 0; j < numberComponents; ++j) {
      vals[j] = data(j) * geoShapes.shapes(i) * dA;
    }
    paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, vals, V.getId(),
                                         numberComponents, name);
  }
  // auto shapes = this->getH1Shapes(pointers, 1, IntegrationPt);
  // std::vector<prec> vals(numberComponents);

  // auto jaco = this->getJacobian(pointers, IntegrationPt);
  // auto dA = jaco.determinant() * IntegrationPt.weight;
  // for (auto i = 0; i < 4; ++i) {
  //   auto &V = pointers.getGeometryData()->getVertex(this->verts[i]);
  //   for (auto j = 0; j < numberComponents; ++j) {
  //     vals[j] = data(j) * shapes.shapes(i) * dA;
  //   }
  //   paraviewAdapter.SumPointDataWeighted(mainMesh, subMesh, vals, V.getId(),
  //                                        numberComponents, name);
  // }
}

void QuadraticQuadrilateral::flip() {
  throw std::runtime_error(
      "ERROR in QuadraticQuadrilateral::flip: Method not implemented!");
}

void QuadraticQuadrilateral::rotate(indexType n) {
  throw std::runtime_error(
      "ERROR in QuadraticQuadrilateral::rotate: Method not implemented!");
}

auto QuadraticQuadrilateral::computeMeanCoordinate(PointerCollection &pointers)
    -> Types::Vector3<prec> {
  throw std::runtime_error("ERROR in QuadraticQuadrilateral::computeMeanCoordinate: "
                           "Method not implemented!");
  return Types::Vector3<prec>();
}

auto QuadraticQuadrilateral::getGeometryShapes(PointerCollection &pointers,
                                               IntegrationPoint &ip)
    -> H1Shapes {
  H1Shapes shapes(9,2);

  auto N1 = [](prec xi) { return (xi - prec(1)) * xi / prec(2); };
  auto N2 = [](prec xi) { return (xi + prec(1)) * xi / prec(2); };
  auto N3 = [](prec xi) { return (prec(1) - xi) * (prec(1) + xi); };
  auto shp = [N1, N2, N3](prec xi) {
    return std::vector<prec>({N1(xi), N2(xi), N3(xi)});
  };

  auto N1xi = [](prec xi) { return xi - prec(0.5); };
  auto N2xi = [](prec xi) { return (xi + prec(0.5)); };
  auto N3xi = [](prec xi) { return -prec(2) * xi; };
  auto shpxi = [N1xi, N2xi, N3xi](prec xi) {
    return std::vector<prec>({N1xi(xi), N2xi(xi), N3xi(xi)});
  };

  auto shapevecXi = shp(ip.xi);
  auto shapevecEta = shp(ip.eta);
  auto DshapevecXi = shpxi(ip.xi);
  auto DshapevecEta = shpxi(ip.eta);

  std::array<indexType, 9> nnXi = {0, 1, 1, 0, 2, 1, 2, 0, 2};
  std::array<indexType, 9> nnEta = {0, 0, 1, 1, 0, 2, 1, 2, 2};

  for (auto i = 0; i < 9; ++i) {
    shapes.shapes(i) = shapevecXi[nnXi[i]] * shapevecEta[nnEta[i]];
    shapes.shapeDeriv(0, i) = DshapevecXi[nnXi[i]] * shapevecEta[nnEta[i]];
    shapes.shapeDeriv(1, i) = shapevecXi[nnXi[i]] * DshapevecEta[nnEta[i]];
  }

  return shapes;
}





const GeometryTypes QuadraticQuadrilateral::type =
    GeometryTypes::QuadraticQuadrilateral;
} // namespace HierAMuS::Geometry
