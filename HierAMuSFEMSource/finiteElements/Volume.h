// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "GenericFiniteElement.h"
#include "geometry/Faces.h"
#include "geometry/Volumes.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

namespace HierAMuS::FiniteElement {
class Volume : public GenericFiniteElement {
  using ptrCol = PointerCollection;

public:
  Volume() = default;
  ~Volume() override;

  auto getType() -> Elementtypes override;

  void setVolume(indexType volumeIn) override;
  auto getVertexIds(PointerCollection& pointers) -> std::vector<indexType> override;

  auto getVertex(ptrCol &pointers, indexType localNumber)
      -> Geometry::Vertex & override;
  auto getEdge(ptrCol &pointers, indexType localNumber)
      -> Geometry::Edges & override;
  auto getFace(ptrCol &pointers, indexType localNumber)
      -> Geometry::Faces * override;
  auto getVolume(ptrCol &pointers, indexType localNumber)
      -> Geometry::Volumes * override;
  

  void setAllNodeBoundaryConditionMeshId(ptrCol &pointers, indexType meshId,
                                         indexType dof) override;

  // Structure
  /**
   * @brief getNumberOfVertices returns the number of vertices.
   * @return indexType, the number of vertices.
   */
  auto getNumberOfVertices(PointerCollection& pointers) -> indexType override;
  /**
   * @brief getNumberOfEdges returns the number of edges.
   * @return indexType, the number of edges.
   */
  auto getNumberOfEdges(PointerCollection& pointers) -> indexType override;
  /**
   * @brief getNumberOfFace returns the number of faces.
   * @return indexType, the number of faces.
   */
  auto getNumberOfFaces(PointerCollection& pointers) -> indexType override;

  // Geometric mapping
  /**
   * @brief Get the Jacobian matrix at the current IntegrationPoint. New version
   *
   * @param pointers, object containtig the pointers to global data.
   * @param IntegrationPt, current IntegrationPoint.
   * @return Types::MatrixXX<prec>, Jacobian matrix at the current IntegrationPoint.
   */
   auto getJacobian(ptrCol &pointers, IntegrationPoint &IntegrationPt)
      -> Types::MatrixXX<prec> override;
  void getJacobian(ptrCol &pointers, Types::Matrix33<prec> &jacobi, prec xsi,
                   prec eta, prec zeta) override;

  // H1 Shapes
  /**
   * @brief Assignes the nodes to the geometry object for H1 shape functions with given order associated with the specific meshId.
   *        A single element can have multiple H1 nodes with different meshids (e.g. for different solutions like displacements and rotations).
   *
   * @param pointers, object containtig the pointers to global data.
   * @param meshId, mesh id of the nodes for the H1 shape functions.
   * @param order, order of the H1 shape functions.
   * @param type, type of the degrees of freedom, currently unused.
   */
  void setH1Shapes(ptrCol &pointers, indexType meshid,
                   indexType order) override;
  void getH1Dofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                 indexType meshID, indexType order) override;
  auto getH1Nodes(ptrCol &pointers, indexType meshID, indexType order)
      -> std::vector<GenericNodes *> override;
  void getH1Shapes(ptrCol &pointers, indexType order,
                           Types::MatrixXX<prec> &jacobi,
                           Types::VectorX<prec> &shape,
                           Types::MatrixXX<prec> &shapeDerivative,
                           IntegrationPoint &IntegrationPt) override;
  void getH1Shapes(ptrCol &pointers, indexType order,
                   Types::Matrix33<prec> &jacobi, Types::VectorX<prec> &shape,
                   Types::Matrix3X<prec> &dshape, prec xi, prec eta, prec zeta) override;

  /**
   * @brief New version of getH1Shapes.
   *
   * @param pointers[in], object containtig the pointers to global data.
   * @param order[in], order of the H1 shape functions.
   * @param IntegrationPt[in], current IntegrationPoint.
   * @return H1Shapes, H1 shape functions evaluation ath the IntegrationPoint.
   */
  auto getH1Shapes(ptrCol &pointers, indexType order,
                           Types::MatrixXX<prec> &jacobi,
                           IntegrationPoint &IntegrationPt)
      -> Geometry::H1Shapes override;

  /**
   * @brief Get the Integration Points object. Returns the IntegrationPoints object, already set to the correct type. Only the integration order needs to be specified.
   *
   * @param[in] pointers, object containtig the pointers to global data.
   * @return IntegrationPoints, IntegrationPoints object set to the correct type.
   */
  auto getIntegrationPoints(ptrCol &pointers) -> IntegrationPoints override;

  // HDiv Shapes
  void setHDivShapes(ptrCol &pointers, indexType meshid,
                     indexType order, NodeTypes type) override;

  void getHDivDofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                   indexType meshID, indexType order) override;
  void getHDivShapes(ptrCol &pointers, indexType order,
                     Types::Matrix22<prec> &jacobi,
                     Types::Matrix2X<prec> &shape, Types::VectorX<prec> &dshape,
                     prec xi, prec eta) override;

  void toParaviewAdapter(PointerCollection &ptrCol, vtkPlotInterface &catalyst,
                         const ParaviewSwitch &ToDo);

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
                            indexType meshId, indexType order, std::string name) override;
  void projectDataToParaviewVertices(
      PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
      indexType mainMesh, indexType subMesh, indexType order,
      IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
      indexType numberComponents, std::string name) override;

protected:
  indexType volume;
};
} // namespace HierAMuS
