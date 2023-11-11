// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <geometry/Faces/FacesDataInterface.h>
#include <types/MatrixTypes.h>
#include <vector>

namespace HierAMuS::Geometry {

class ScaledBoundary2DData
    : public FacesDataInterface<-1, -1, ScaledBoundary2DData> {
public:
  ScaledBoundary2DData();
  ~ScaledBoundary2DData() override;

  auto getRuntimeObject(GeometryData &geoData)
      -> std::shared_ptr<FacesRuntime> override;

  auto getType() -> const GeometryTypes & override;
  auto getNumberOfVerts() -> indexType override {
    return this->m_verts.size();
  };
  /**
   * @brief Returns the number of edges of the element.
   * @return Number of edges as an indexType
   */
  auto getNumberOfEdges() -> indexType override {
    return this->m_edges.size();
  };
  void setVerts(GeometryData &geoData,
                std::vector<indexType> &vertsIn) override;
  void setEdges(const std::vector<indexType> &edgesIn) override;
  // void setFaces(std::vector<indexType>& facesIn);
  auto getVertexNumbers() -> std::vector<indexType> override;
  void getVerts(std::vector<VertexData *> &vertsOut) override;
  void setScalingCenter(prec x, prec y);
  void computeScalingCenter();
  auto getEdgeNumbers() -> std::vector<indexType> override;
  void getEdgeNumbers(std::vector<indexType> &edgesOut) override {
    edgesOut = std::vector<indexType>(m_edges.begin(), m_edges.end());
  };
  void getEdges(std::vector<EdgesData *> &edgesOut) override;
  auto getEdgeNumber(indexType num) -> indexType override {
    return this->m_edges[num];
  };
  /** @brief Returns the global edge number based on the local edge number of
   * the face.
   *
   *  detailed description
   *
   *  @param[in] local_number local edge number of the face
   *  @return indexType global edge number
   */
  auto getEdge(indexType local_number) -> EdgesData * override;

  /** @brief Check if the face has the three vertices.
   *
   * @param[in] v1 Global number of a corner vertex.
   * @param[in] v2 Global number of a corner vertex.
   * @param[in] v3 Global number of a corner vertex.
   * @return true If the face has the three vertices.
   * @return false If the face does not have the three vertices.
   */
  auto hasVertices(indexType v1, indexType v2, indexType v3) -> bool override;

  void print(spdlog::logger &Logger) override;

  auto getCoordinates(prec xi, prec eta)
      -> Types::Vector3<prec> override;

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
  void getH1DofsInternal(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) override;

  auto getH1NodesList(indexType meshID,
                      indexType order) -> MeshIdNodeList override;

  void getH1Shapes(indexType order,
                   Types::VectorX<prec> &shape,
                   Types::Matrix2X<prec> &shapeDerivative, prec xsi,
                   prec eta) override;

  void getH1ShapesInternal(indexType order,
                           Types::VectorX<prec> &shape,
                           Types::Matrix2X<prec> &shapeDerivative, prec xsi,
                           prec eta) override;

  auto getH1Shapes(indexType order,
                   IntegrationPoint &IntegrationPt) -> H1Shapes override;

  // HDivShapes
  void setHDivShapes(indexType meshid,
                     indexType order, NodeTypes type) override;

  void getHDivDofs(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                   indexType order, NodeTypes type) override;
  void getHDivShapes(indexType order,
                     Types::Matrix2X<prec> &shape, Types::VectorX<prec> &dshape,
                     prec xi, prec eta) override;

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

  static auto getName() -> std::string {
    return "Scaled Boundary Face Element";
  };

private:
  // indexType scalingCenter;
  std::vector<prec> scalingcoor;
  std::vector<indexType> m_edges;
  std::vector<indexType> m_verts;

  std::vector<EdgesData *> m_edges_pointers;
  std::vector<VertexData *> m_vert_pointers;
  static const GeometryTypes type;
};
} // namespace HierAMuS::Geometry
