// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause




namespace HierAMuS::Geometry {
	enum class GeometryTypes{Generic,
		Vertex,
		Edges=100,LinearEdge=102,QuadraticEdge=103,
		Faces=200,LinearTriangle=203,LinearQuadrilateral=204, QuadrilateralNodal,
		Volumes=300,Tetra=303
	};

}
