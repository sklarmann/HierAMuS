// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <array>

#include "geometry/Volumes/VolumesH1Interface.h"
#include <geometry/GeometryTypes.h>
#include <geometry/Volumes/VolumesRuntimeDataInterface.h>

#include <types/MatrixTypes.h>

#include <vector>

namespace HierAMuS::Geometry {
class LinearPrismData;
class LinearPrismRuntime
    : public VolumesRuntimeDataInterface<6, 9, 5, LinearPrismData,
                                         LinearPrismRuntime>,
      public VolumesH1Interface {
public:
  LinearPrismRuntime(GeometryData &geoData, LinearPrismData &data_element);
  ~LinearPrismRuntime();

  auto getH1Volume() -> VolumesH1Interface * override { return this; };

  auto getType() -> const GeometryTypes & override;
  void print(spdlog::logger &Log) override;

  void getEdgeNumbers(std::vector<indexType> &edgesOut) override;
  void getFaceNumbers(std::vector<indexType> &facesOut) override;

  virtual auto getCoordinates(IntegrationPoint &IntPoint)
      -> Types::Vector3<prec> override;

  // Geometric mapping

  // H1 Shapes
  void setH1Shapes(indexType meshId,
                   indexType order, NodeTypes type) override;
  void setH1ShapesInternal(indexType meshId,
                           indexType order, NodeTypes type) override;

  void getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                 indexType order) override;

  void getH1DofsInternal(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) override;

  void getH1Shapes(indexType order,
                   Types::VectorX<prec> &shape,
                   Types::Matrix3X<prec> &shapeDerivative, prec xsi, prec eta,
                   prec zeta) override;

  void getH1ShapesInternal(indexType order,
                           Types::VectorX<prec> &shape,
                           Types::Matrix3X<prec> &shapeDerivative, prec xsi,
                           prec eta, prec zeta) override;

  
  
  /**
   * @brief New version of getH1Shapes.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param order[in], order of the H1 shape functions.
   * @param IntegrationPt[in], integration point.
   * @return H1Shapes, H1 shape functions.
   */
  auto getH1Shapes(indexType order,
                   IntegrationPoint &IntegrationPt) -> H1Shapes override;
  /**
   * @brief New version of getH1ShapesInternal.
   * Computes and returns the H1 bubble functions of the volume element.
   *
   * @todo Implementation of higher order H1 shapes.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param order[in], order of the H1 shape functions.
   * @param IntegrationPt[in], integration point.
   * @return H1Shapes, H1 shape functions.
   */
  auto getH1ShapesInternal(indexType order,
                           IntegrationPoint &IntegrationPt)
      -> H1Shapes override;

  
  auto getH1Nodes(indexType meshID, indexType order)
      -> std::vector<GenericNodes *> override;
  auto getH1NodesInternal(indexType meshID, indexType order)
      -> std::vector<GenericNodes *> override;

private:
  const static indexType m_numberOfVerts = 6;
  const static indexType m_numberOfEdges = 9;
  const static indexType m_numberOfFaces = 5;

  LinearPrismData &m_LinearPrism_data;
  std::array<VertexRuntime *, m_numberOfVerts> m_Vertices;
  std::array<std::shared_ptr<EdgesRuntime>, m_numberOfEdges> m_Edges;
  std::array<std::shared_ptr<FacesRuntime>, m_numberOfFaces> m_Faces;

  std::array<indexType, 9> edges;
  std::array<indexType, 5> faces;

  // indexType verts[6], edges[9], faces[5];
  static const GeometryTypes type;
};

} // namespace HierAMuS::Geometry