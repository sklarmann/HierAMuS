// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "equations/DegreeOfFreedom.h"
#include "equations/GenericNodes.h"
#include "pointercollection/pointercollection.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <forwarddeclaration.h>
#include <geometry/GeometryData.h>
#include <geometry/Special.h>

#include <geometry/GeometryTypes.h>
#include <types/MatrixTypes.h>

#include <map>
#include <vector>

namespace HierAMuS::Geometry {

class BeamInterface3D : public Special {
public:
  BeamInterface3D();
  ~BeamInterface3D() override;
  auto getGroupType() -> const GeometryTypes & override { return HierAMuS::Geometry::BeamInterface3D::type; };

  
  /**
   * @brief Returns the number of edges of the element.
   * @return Number of edges as an indexType
   */
  auto getNumberOfEdges() -> indexType override { return 0; };

  /**
   * @brief Sets the faces of the interface. There can be an arbitrary amount of
   * faces, depending on the cross-section of the beam.
   *
   * @param faces[in], vector of global ids of the faces.
   */
  void setFaces(std::vector<indexType> &faces) override;

  void setBeamVertex(indexType vertex) override;

  void setH1Shapes(PointerCollection &pointers, indexType meshId,
                           indexType order, NodeTypes type) override;
  void setH1ShapesBeamRot(PointerCollection &pointers, indexType meshId, NodeTypes type);
  /**
   * @brief Get the Integration Points for the geometry element. Currently the
   * IntegrationType::Scaled2D scheme is used. Number of faces are treated as
   * number of sections.
   *
   * @param pointers[in], pointers to global data.
   * @return IntegrationPoints, integration points for the geometry element.
   */
  auto getIntegrationPoints(PointerCollection &pointers, indexType elementId)
  -> IntegrationPoints override;

  void computeWarpingShapes(PointerCollection &pointers);

  auto getCoordinates(PointerCollection& pointers, IntegrationPoint& IntPoint) -> Types::Vector3<prec> override;

  
  void print(PointerCollection &pointers) override {};

private:
  void createNodeShapeMapping(PointerCollection &pointers);
  void computeGeometry(PointerCollection &pointers);
  auto getLocalCoordinate(PointerCollection &pointers,
                          IntegrationPoint &integrationPoint) -> Types::Vector3<prec>;
  auto static getFaceIntegrationPoint(IntegrationPoint &integrationPoint) -> IntegrationPoint;


  static const GeometryTypes type;

  // Warping shape data
  std::map<indexType, indexType> m_nodeShapeMapping;
  indexType m_H1MeshId; // Displacement nodes will be used to sort the warping shapes
  indexType m_warpOrder;
  indexType m_numberOfWarpingShapes;

  // Geometry data
  prec m_length;
  Types::Vector3<prec> m_projectedCoordinate;
  Types::Vector3<prec> m_A1;
  Types::Vector3<prec> m_A2;
  Types::Vector3<prec> m_A3;

  // Description
  std::vector<indexType> m_faces;
  indexType m_beamVertex;

  prec A;

};

} // namespace HierAMuS
