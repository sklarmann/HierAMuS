// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause




#pragma once



namespace HierAMuS::Geometry {
enum class [[nodiscard]] GeometryTypes {
  Generic,
  Vertex,
  Edges = 100,
  LinearEdge = 102,
  QuadraticEdge = 103,
  Faces = 200,
  LinearTriangle = 203,
  LinearQuadrilateral = 204,
  QuadraticQuadrilateral = 205,
  QuadrilateralNodal,
  Volumes = 300,
  Tetra = 303,
  LinearBrick = 308,
  LinearPrism = 306,
  Special = 400,
  BeamInterface2D = 401,
  BeamInterface3D = 402,
  ScaledBoundary2D = 502,
};

enum class [[nodiscard]] ShapeFunctionTypes { H0, H1, HDiv, HCurl };

} // namespace HierAMuS