// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <geometry/Faces/FacesRuntime.h>
#include "geometry/Faces/FacesH1Interface.h"
#include <types/MatrixTypes.h>
#include <vector>

namespace HierAMuS::Geometry {
class ScaledBoundary2DData;

class ScaledBoundary2DRuntime : public FacesRuntime, public FacesH1Interface {
public:
  ScaledBoundary2DRuntime(GeometryData &geoData,
                          ScaledBoundary2DData &base_element);
  ~ScaledBoundary2DRuntime() override;

  auto getH1Face() -> FacesH1Interface * override { return this; };

  auto getType() -> const GeometryTypes & override;
  auto getNumberOfVerts() -> indexType override { return this->verts.size(); };
  /**
   * @brief Returns the number of edges of the element.
   * @return Number of edges as an indexType
   */
  auto getNumberOfEdges() -> indexType override { return this->edges.size(); };
  // void setFaces(std::vector<indexType>& facesIn);
  auto getVertexNumbers() -> std::vector<indexType> override;
  void getVertices(std::vector<VertexRuntime *> &vertsOut) override;
  auto getVertex(indexType local_number) -> Geometry::VertexRuntime * override;
  void setScalingCenter(prec x, prec y);
  void computeScalingCenter();
  auto getEdgeNumbers() -> std::vector<indexType> override;
  void getEdges(std::vector<EdgesRuntime *> &edgesOut) override;

  /** @brief Returns the global edge number based on the local edge number of
   * the face.
   *
   *  detailed description
   *
   *  @param[in] local_number local edge number of the face
   *  @return indexType global edge number
   */
  auto getEdgeNumber(indexType local_number) -> indexType override;
  auto getEdge(indexType local_number) -> Geometry::EdgesRuntime * override;

  /** @brief Check if the face has the three vertices.
   *
   * @param[in] v1 Global number of a corner vertex.
   * @param[in] v2 Global number of a corner vertex.
   * @param[in] v3 Global number of a corner vertex.
   * @return true If the face has the three vertices.
   * @return false If the face does not have the three vertices.
   */
  auto hasVertices(indexType v1, indexType v2, indexType v3) -> bool override;

  void print(spdlog::logger &Log) override;

  auto getCoordinates(IntegrationPoint &integrationPoint)
      -> Types::Vector3<prec> override;

  auto getIntegrationPoints(indexType elementId) -> IntegrationPoints override;

  auto getJacobian(IntegrationPoint &IntegrationPt)
      -> Types::Matrix22<prec> override;

  // H1Shapes
  void setH1Shapes(indexType meshId,
                   indexType order, NodeTypes type) override;
  void setH1ShapesInternal(indexType meshId,
                           indexType order, NodeTypes type) override;
  void getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                 indexType order) override;
  auto getH1Dofs(indexType meshID, indexType order)
      -> std::vector<DegreeOfFreedom *> override;
  void getH1DofsInternal(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) override;

  auto getH1NodesList(indexType meshID,
                      indexType order) -> MeshIdNodeList override;

  auto getH1Shapes(indexType order, IntegrationPoint &IntegrationPt)
      -> H1Shapes override;
  auto getH1ShapesInternal(indexType order, IntegrationPoint &IntegrationPt,
                           faceorientation orientation = faceorientation::p_1)
      -> H1Shapes override;
  auto getH1Nodes(indexType meshID, indexType order)
      -> std::vector<GenericNodes *> override;
  auto getH1NodesInternal(indexType meshID, indexType order)
      -> std::vector<GenericNodes *> override;



  void setAllNodeBoundaryConditionMeshId(indexType meshId,
                                         indexType dof) override;

  // L2 shapes

  void getL2Dofs(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                 indexType order) override;

  auto getL2Shapes(indexType order,
                   IntegrationPoint &IntegrationPt) -> L2Shapes override;

  void getL2ShapesInternal(indexType order,
                           Types::VectorX<prec> &shape,
                           Types::Matrix2X<prec> &shapeDerivative, prec xsi,
                           prec eta) override;

  void setL2Shapes(indexType meshId,
                   indexType order, NodeTypes type) override;

  void geometryToParaview(vtkPlotInterface &paraviewAdapter, indexType mainMesh,
                          indexType subMesh) override;

  void computeWeightsParaview(vtkPlotInterface &paraviewAdapter,
                              indexType mainMesh, indexType subMesh) override;

  void H1SolutionToParaview(vtkPlotInterface &paraviewAdapter,
                            indexType mainMesh, indexType subMesh,
                            indexType order, Types::VectorX<prec> &solution,
                            std::string &name) override;

  void H1DataToParaview(vtkPlotInterface &paraviewAdapter, indexType mainMesh,
                        indexType subMesh, Types::VectorX<prec> &Data,
                        indexType numberComponents, indexType order,
                        std::string &name) override;

  void projectDataToParaviewVertices(
      vtkPlotInterface &paraviewAdapter,
      indexType mainMesh, indexType subMesh, indexType order,
      IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
      indexType numberComponents, std::string name) override;

  void checkUpdateElement(EquationHandler &eqHandler,
                          GeometryData &geoData) override{};

  /**
   * @brief Changes orientation of the face (reverts the face normal vector).
   */
  void flip() override;

  /**
   * @brief Rotates the face n times.
   *
   * @param [in] n Number of rotations.
   */
  void rotate(indexType n) override;
  /**
   * @brief Computes the mean coordinate of the face element.
   *
   * @param [in] pointers Pointer collection to global data.
   *
   * @return Mean coordinate of the face element.
   */
  auto computeMeanCoordinate()
      -> Types::Vector3<prec> override;

  /**
   * @brief Get the Face Normal vector
   * Only virtual method, must be implemented in the derived face classes.
   *
   * @param pointers[in], pointers to global data.
   * @return Types::Vector3<prec>, Normal vector of the face
   */
  auto getFaceNormal()
      -> Types::Vector3<prec> override;

  void set_geometry_pointers(GeometryData &geoData) override;

  auto getVertexNumber(indexType localNumber) -> indexType override;

private:
  ScaledBoundary2DData &m_Scaled_data;

  // indexType scalingCenter;
  std::vector<prec> scalingcoor;
  std::vector<indexType> edges;
  std::vector<indexType> verts;

  std::vector<EdgesData *> m_edge_pointers;
  std::vector<VertexData *> m_vert_pointers;

  std::vector<VertexRuntime *> m_Verts;
  std::vector<std::shared_ptr<EdgesRuntime>> m_Edges;
  static const GeometryTypes type;
};
} // namespace HierAMuS::Geometry
