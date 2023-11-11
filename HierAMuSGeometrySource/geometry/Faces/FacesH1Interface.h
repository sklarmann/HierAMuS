// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "datatypes.h"
#include "geometry/GeometryShape.h"
#include "geometry/Faces/FaceOrientationFlags.h"
#include "MeshIdNodeList.h"
#include "Nodetypes.h"

namespace HierAMuS {
class GenericNodes;
class DegreeOfFreedom;
} // namespace HierAMuS

namespace HierAMuS::Geometry {

class FacesH1Interface {
public:
  // FacesH1Interface();
  ~FacesH1Interface() = default;

  virtual auto getH1Face() -> FacesH1Interface * = 0;

  virtual void setH1Shapes(indexType meshId, indexType order,
                           NodeTypes type) = 0;
  virtual void setH1ShapesInternal(indexType meshId, indexType order,
                                   NodeTypes type) = 0;

  virtual auto getH1Dofs(indexType meshID, indexType order)
      -> std::vector<DegreeOfFreedom *> = 0;

  virtual void getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) = 0;

  virtual void getH1DofsInternal(std::vector<DegreeOfFreedom *> &Dofs,
                                 indexType meshID, indexType order) = 0;

  virtual auto getH1NodesList(indexType meshID, indexType order)
      -> MeshIdNodeList = 0;

  virtual auto getH1Nodes(indexType meshID, indexType order)
      -> std::vector<GenericNodes *> = 0;

  virtual auto getH1NodesInternal(indexType meshID, indexType order)
      -> std::vector<GenericNodes *> = 0;

  virtual auto getH1Shapes(indexType order, IntegrationPoint &IntegrationPt)
      -> H1Shapes = 0;

  virtual auto
  getH1ShapesInternal(indexType order, IntegrationPoint &IntegrationPt,
                      faceorientation orientation = faceorientation::p_1)
      -> H1Shapes = 0;
};

} // namespace HierAMuS::Geometry
