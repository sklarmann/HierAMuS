// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <forwarddeclaration.h>

#include <geometry/Faces.h>
#include <types/MatrixTypes.h>
#include <vector>

namespace HierAMuS::Geometry {

class ScaledBoundary2D : public Faces {
public:
  ScaledBoundary2D();
  ~ScaledBoundary2D() override;
  auto getType() -> const GeometryTypes & override;
  auto getNumberOfVerts() -> indexType override { return this->verts.size(); };
  /**
   * @brief Returns the number of edges of the element.
   * @return Number of edges as an indexType
   */
  auto getNumberOfEdges() -> indexType override { return this->edges.size(); };
  void setVerts(GeometryData &geoData, std::vector<indexType> &vertsIn) override;
  void setEdges(const std::vector<indexType> &edgesIn) override;
  // void setFaces(std::vector<indexType>& facesIn);
  void getVerts(std::vector<indexType> &vertsOut) override;
  void getVerts(PointerCollection &pointers,
                std::vector<Base *> &vertsOut) override;
  void setScalingCenter(prec x, prec y);
  void computeScalingCenter(PointerCollection &pointers);
  void getEdges(std::vector<indexType> &edgesOut) override;
  void getEdges(PointerCollection &pointers,
                std::vector<Base *> &edgesOut) override;

  /** @brief Returns the global edge number based on the local edge number of the face.
   *
   *  detailed description
   *
   *  @param[in] local_number local edge number of the face
   *  @return indexType global edge number
   */
  auto getEdge(indexType local_number) -> indexType override;

  /** @brief Check if the face has the three vertices.
    *
    * @param[in] v1 Global number of a corner vertex.
    * @param[in] v2 Global number of a corner vertex.
    * @param[in] v3 Global number of a corner vertex.
    * @return true If the face has the three vertices.
    * @return false If the face does not have the three vertices.
    */
  auto hasVertices(indexType v1, indexType v2, indexType v3) -> bool override;

  void print(PointerCollection &pointers) override;

  auto getCoordinates(PointerCollection &pointers, prec xi,
                                      prec eta) -> Types::Vector3<prec> override;

  auto getIntegrationPoints(PointerCollection &pointers, indexType elementId)
  -> IntegrationPoints override;

  auto getJacobian(PointerCollection &pointers, IntegrationPoint &IntegrationPt)
      -> Types::MatrixXX<prec> override;

  // H1Shapes
  void setH1Shapes(PointerCollection &pointers, indexType meshId,
                   indexType order, NodeTypes type) override;
  void setH1ShapesInternal(PointerCollection &pointers, indexType meshId,
                           indexType order, NodeTypes type) override;
  void getH1Dofs(PointerCollection &pointers,
                 std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                 indexType order) override;
  void getH1DofsInternal(PointerCollection &pointers,
                         std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) override;

  void getH1Shapes(PointerCollection &pointers, indexType order,
                   Types::VectorX<prec> &shape,
                   Types::Matrix2X<prec> &shapeDerivative, prec xsi,
                   prec eta) override;

  void getH1ShapesInternal(PointerCollection &pointers, indexType order,
                           Types::VectorX<prec> &shape,
                           Types::Matrix2X<prec> &shapeDerivative, prec xsi,
                           prec eta) override;

  auto getH1Shapes(PointerCollection &pointers, indexType order,
                   IntegrationPoint &IntegrationPt) -> H1Shapes override;

  // HDivShapes
  void setHDivShapes(PointerCollection &pointers, indexType meshid,
                     indexType order, NodeTypes type) override;

  void getHDivDofs(PointerCollection &pointers,
                   std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                   indexType order, NodeTypes type) override;
  void getHDivShapes(PointerCollection &pointers, indexType order,
                     Types::Matrix2X<prec> &shape, Types::VectorX<prec> &dshape,
                     prec xi, prec eta) override;

   void setAllNodeBoundaryConditionMeshId(PointerCollection &pointers,
                                                 indexType meshId,
                                                 indexType dof) override;

  // L2 shapes

  void getL2Dofs(PointerCollection &pointers,
                 std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                 indexType order) override;

  auto getL2Shapes(PointerCollection &pointers, indexType order,
                   IntegrationPoint &IntegrationPt) -> L2Shapes override;

  void getL2ShapesInternal(PointerCollection &pointers, indexType order,
                           Types::VectorX<prec> &shape,
                           Types::Matrix2X<prec> &shapeDerivative, prec xsi,
                           prec eta) override;

  void setL2Shapes(PointerCollection &pointers, indexType meshId,
                   indexType order, NodeTypes type) override;

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
                            std::string &name) override;

  void H1DataToParaview(PointerCollection &pointers,
                        vtkPlotInterface &paraviewAdapter, indexType mainMesh,
                        indexType subMesh, Types::VectorX<prec> &Data,
                        indexType numberComponents, indexType order,
                        std::string &name) override;

  void projectDataToParaviewVertices(
      PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
      indexType mainMesh, indexType subMesh, indexType order,
      IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
      indexType numberComponents, std::string name) override;

  void checkUpdateElement(GeometryData &geoData) override{};


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
  auto computeMeanCoordinate(PointerCollection &pointers)
      -> Types::Vector3<prec> override;

private:
  // indexType scalingCenter;
  std::vector<prec> scalingcoor;
  std::vector<indexType> edges;
  std::vector<indexType> verts;
  static const GeometryTypes type;
};
} // namespace HierAMuS
