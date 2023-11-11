// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "GenericFiniteElementInterface.h"
#include "geometry/Faces/FacesData.h"
#include "geometry/Volumes/VolumesData.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "geometry/VertexRuntime.h"

namespace HierAMuS::Geometry{
class VolumesRuntime;
}

namespace HierAMuS::FiniteElement {
class VolumeConstraint
    : public GenericFiniteElementInterface<VolumeConstraint> {
  using ptrCol = PointerCollection;

public:
  VolumeConstraint() = default;
  ~VolumeConstraint() override;
  void set_pointers(PointerCollection &pointers) override;

  void setVerts(std::vector<indexType> &vertsIn) override;
  void setVolume(indexType volumeIn) override;

  auto getVertexCoordinates(ptrCol &pointers) -> Types::Vector3<prec>;
  auto getVolumeCoordinates(ptrCol &pointers, IntegrationPoint &IntegrationPt)
      -> Types::Vector3<prec>;
  
  void setVertexNodes(PointerCollection &pointers, indexType meshId);
  void getVertexDofs(PointerCollection &pointers,
                     std::vector<DegreeOfFreedom *> &Dofs, indexType meshID);

  auto getVertexNodes(PointerCollection &pointers, indexType meshId)
      -> std::vector<GenericNodes *>;

  // Geometric mapping
  /**
   * @brief Get the Jacobian matrix at the current IntegrationPoint. New version
   *
   * @param pointers, object containtig the pointers to global data.
   * @param IntegrationPt, current IntegrationPoint.
   * @return Types::MatrixXX<prec>, Jacobian matrix at the current
   * IntegrationPoint.
   */
  auto getJacobian(ptrCol &pointers, IntegrationPoint &IntegrationPt)
      -> Types::Matrix33<prec>;
  // H1 Shapes
  /**
   * @brief Assignes the nodes to the geometry object for H1 shape functions
   * with given order associated with the specific meshId. A single element can
   * have multiple H1 nodes with different meshids (e.g. for different solutions
   * like displacements and rotations).
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
  /**
   * @brief New version of getH1Shapes.
   *
   * @param pointers[in], object containtig the pointers to global data.
   * @param order[in], order of the H1 shape functions.
   * @param IntegrationPt[in], current IntegrationPoint.
   * @return H1Shapes, H1 shape functions evaluation ath the IntegrationPoint.
   */
  auto getH1Shapes(ptrCol &pointers, indexType order,
                   Types::Matrix33<prec> &jacobi,
                   IntegrationPoint &IntegrationPt) -> Geometry::H1Shapes;
  /**
   * @brief Get the Integration Points object. Returns the IntegrationPoints
   * object, already set to the correct type. Only the integration order needs
   * to be specified.
   *
   * @param[in] pointers, object containtig the pointers to global data.
   * @return IntegrationPoints, IntegrationPoints object set to the correct
   * type.
   */
  auto getIntegrationPoints(ptrCol &pointers) -> IntegrationPoints override;

private:
  indexType m_sharedVertex;
  indexType m_volume;

  Geometry::VolumesData *m_volume_pointer;
  Geometry::VertexData *m_vertex_pointer;

  Geometry::VertexRuntime* m_Vertex;
  std::shared_ptr<Geometry::VolumesRuntime> m_Volume;
};
} // namespace HierAMuS
