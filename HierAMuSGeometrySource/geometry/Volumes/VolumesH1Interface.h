// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "datatypes.h"
#include "geometry/GeometryShape.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "MeshIdNodeList.h"
#include "Nodetypes.h"

namespace HierAMuS {
class GenericNodes;
class DegreeOfFreedom;
} // namespace HierAMuS

namespace HierAMuS::Geometry {

class VolumesH1Interface {
public:
  // VolumesH1Interface();
  ~VolumesH1Interface() = default;

  virtual auto getH1Volume() -> VolumesH1Interface * = 0;

  virtual void setH1Shapes(indexType meshId, indexType order,
                           NodeTypes type) = 0;
  virtual void setH1ShapesInternal(indexType meshId, indexType order,
                                   NodeTypes type) = 0;

  virtual void getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) = 0;

  virtual void getH1DofsInternal(std::vector<DegreeOfFreedom *> &Dofs,
                                 indexType meshID, indexType order) = 0;

  virtual void getH1Shapes(indexType order, Types::VectorX<prec> &shape,
                           Types::Matrix3X<prec> &shapeDerivative, prec xsi,
                           prec eta, prec zeta) = 0;

  virtual void getH1ShapesInternal(indexType order, Types::VectorX<prec> &shape,
                                   Types::Matrix3X<prec> &shapeDerivative,
                                   prec xsi, prec eta, prec zeta) = 0;

  /**
   * @brief New version of getH1Shapes.
   * Virtual method to get the H1 shapes of the element.
   *
   * @param order[in], order of the H1 shape functions.
   * @param IntegrationPt[in], integration point.
   * @return H1Shapes, H1 shape functions.
   */
  virtual auto getH1Shapes(indexType order,
                           IntegrationPoint &IntegrationPt) -> H1Shapes = 0;
  /**
   * @brief New version of getH1ShapesInternal.
   * Computes and returns the H1 bubble functions of the volume element.
   *
   * @param order[in], order of the H1 shape functions.
   * @param IntegrationPt[in], integration point.
   * @return H1Shapes, H1 shape functions.
   */
  virtual auto getH1ShapesInternal(indexType order,
                                   IntegrationPoint &IntegrationPt)
      -> H1Shapes = 0;

  /**
   * @brief Returns a vector of the H1 nodes for the given meshId and shape
   * function order.
   *
   * Here it is a virtual method overridden by the actual geometry object.
   *
   * @param meshID, mesh id of the nodes for the H1 shape functions.
   * @param order, order of the H1 shape functions.
   * @return std::vector<GenericNodes *>, vector of the GenericNodes.
   */
  virtual auto getH1Nodes(indexType meshID, indexType order)
      -> std::vector<GenericNodes *> = 0;
  virtual auto getH1NodesInternal(indexType meshID, indexType order)
      -> std::vector<GenericNodes *> = 0;
};

} // namespace HierAMuS::Geometry
