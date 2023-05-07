// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once

#include "Face.h"
#include "MatrixTypes.h"
#include "pointercollection/pointercollection.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

namespace HierAMuS::FiniteElement {
class FaceConstraint : public Face {
  using ptrCol = PointerCollection;

public:
  FaceConstraint() = default;
  ~FaceConstraint() override;

  auto getType() -> Elementtypes override;

  void setVerts(std::vector<indexType> &vertsIn) override;
  void setFace(indexType faceIn) override;
  auto getVertexIds(PointerCollection &pointers)
      -> std::vector<indexType> override;

  auto getVertex(ptrCol &pointers, indexType localNumber)
      -> Geometry::Vertex & override;
  auto getEdge(ptrCol &pointers, indexType localNumber)
      -> Geometry::Edges & override;
  auto getFace(ptrCol &pointers, indexType localNumber)
      -> Geometry::Faces * override;

  auto getLocalCoordinateSystem(PointerCollection &pointers)
      -> Types::Matrix33<prec>;

  void setAllNodeBoundaryConditionMeshId(ptrCol &pointers, indexType meshId,
                                         indexType dof) override;

  // Structure
  auto getNumberOfVertices(PointerCollection &pointers) -> indexType override;
  auto getNumberOfEdges(PointerCollection &pointers) -> indexType override;

  // Geometric mapping
  auto getJacobian(ptrCol &pointers, IntegrationPoint &IntegrationPt)
      -> Types::MatrixXX<prec> override;
  void getJacobian(ptrCol &pointers, Types::Matrix22<prec> &jacobi, prec xsi,
                   prec eta) override;

  // H1 Shapes
  void setH1Shapes(ptrCol &pointers, indexType meshid,
                   indexType order) override;
  void getH1Dofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                 indexType meshID, indexType order) override;
  auto getH1Nodes(ptrCol &pointers, indexType meshID, indexType order)
      -> std::vector<GenericNodes *> override;
  void getH1Shapes(ptrCol &pointers, indexType order,
                   Types::MatrixXX<prec> &jacobi, Types::VectorX<prec> &shape,
                   Types::MatrixXX<prec> &shapeDerivative,
                   IntegrationPoint &IntegrationPt) override;
  void getH1Shapes(ptrCol &pointers, indexType order,
                   Types::Matrix22<prec> &jacobi, Types::VectorX<prec> &shape,
                   Types::Matrix2X<prec> &dshape, prec xi, prec eta) override;

  auto getH1Shapes(ptrCol &pointers, indexType order,
                   Types::MatrixXX<prec> &jacobi,
                   IntegrationPoint &IntegrationPt)
      -> Geometry::H1Shapes override;

  auto getIntegrationPoints(ptrCol &pointers) -> IntegrationPoints override;

  // Paraview
  void geometryToParaview(PointerCollection &pointers,
                          vtkPlotInterface &paraviewAdapter, indexType mainMesh,
                          indexType subMesh) override;
  void computeWeightsParaview(PointerCollection &pointers,
                              vtkPlotInterface &paraviewAdapter,
                              indexType mainMesh, indexType subMesh) override;
  void H1SolutionToParaview(PointerCollection &pointers,
                            vtkPlotInterface &paraviewAdapter,
                            indexType mainMesh, indexType subMesh,
                            indexType meshId, indexType order,
                            std::string name) override;
  void projectDataToParaviewVertices(
      PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
      indexType mainMesh, indexType subMesh, indexType order,
      IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
      indexType numberComponents, std::string name) override;

  auto getVertexCoordinates(ptrCol &pointers) -> Types::Vector3<prec>;
  auto getFaceCoordinates(ptrCol &pointers, IntegrationPoint &IntegrationPt)
      -> Types::Vector3<prec>;

  void setVertexNodes(PointerCollection &pointers, indexType meshId);
  void getVertexDofs(PointerCollection &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                     indexType meshID);

  auto getTangentG1(PointerCollection &pointers, IntegrationPoint &IntegrationPt)
      -> Types::Vector3<prec>;
  auto getTangentG2(PointerCollection &pointers, IntegrationPoint &IntegrationPt)
      -> Types::Vector3<prec>;

private:
  indexType vertex;
};
} // namespace HierAMuS::FiniteElement
