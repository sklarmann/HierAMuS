// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <array>
#include <forwarddeclaration.h>
#include <geometry/GeometryTypes.h>
#include <geometry/Volumes.h>

#include <types/MatrixTypes.h>

#include <vector>

namespace HierAMuS::Geometry {

class LinearPrism : public Volumes {
public:
  LinearPrism();
  ~LinearPrism();
  auto getType() -> const GeometryTypes & override;
  void print(PointerCollection &pointers) override;

  void setVerts(GeometryData &geoData, std::vector<indexType> &vertsIn) override;
  void setEdges(const std::vector<indexType> &edgesIn) override;
  void setFaces(std::vector<indexType> &facesIn) override;

  void getVerts(std::vector<indexType> &vertsOut) override;
  void getVerts(PointerCollection &pointers, std::vector<Base *> &vertsOut) override;
  auto getNumberOfVerts() -> indexType override { return 8; };
  void getEdges(std::vector<indexType> &edgesOut) override;
  void getEdges(PointerCollection &pointers, std::vector<Base *> &edgesOut) override;
  auto getNumberOfEdges() -> indexType override { return 12; };
  void getFaces(std::vector<indexType> &facesOut) override;
  void getFaces(PointerCollection &pointers, std::vector<Base *> &facesOut) override;

  // Geometric mapping
  void getJacobian(PointerCollection &pointers, Types::Matrix33<prec> &jacobi,
                   prec xsi, prec eta, prec zeta) override;

  // H1 Shapes
  void setH1Shapes(PointerCollection &pointers, indexType meshId,
                   indexType order, NodeTypes type) override;
  void setH1ShapesInternal(PointerCollection &pointers, indexType meshId,
                           indexType order, NodeTypes type) override;

   void getH1Dofs(PointerCollection &pointers,
                         std::vector<DegreeOfFreedom *> &Dofs,
                         indexType meshID, indexType order) override;

  void getH1DofsInternal(PointerCollection &pointers,
                         std::vector<DegreeOfFreedom *> &Dofs,
                         indexType meshID, indexType order) override;

  void getH1Shapes(PointerCollection &pointers, indexType order,
                   Types::VectorX<prec> &shape,
                   Types::Matrix3X<prec> &shapeDerivative, prec xsi,
                   prec eta, prec zeta) override;

  void getH1ShapesInternal(PointerCollection &pointers, indexType order,
                           Types::VectorX<prec> &shape,
                           Types::Matrix3X<prec> &shapeDerivative,
                           prec xsi, prec eta, prec zeta) override;

private:
  std::array<indexType, 6> verts;
  std::array<indexType, 9> edges;
  std::array<indexType, 5> faces;
  
  //indexType verts[6], edges[9], faces[5];
  static const GeometryTypes type;
};

} // namespace HierAMuS