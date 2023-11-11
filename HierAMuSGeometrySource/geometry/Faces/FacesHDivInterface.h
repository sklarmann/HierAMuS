// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "datatypes.h"
#include "geometry/GeometryShape.h"
#include "geometry/Faces/FaceOrientationFlags.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "MeshIdNodeList.h"
#include "Nodetypes.h"

namespace HierAMuS {
class GenericNodes;
class DegreeOfFreedom;
} // namespace HierAMuS

namespace HierAMuS::Geometry {

class FacesHDivInterface {
public:
  // FacesHDivInterface();
  ~FacesHDivInterface() = default;

  virtual auto getHDivFace() -> FacesHDivInterface * = 0;

  

  virtual auto getHDivNodes(indexType meshID,
                            indexType order) -> std::vector<GenericNodes *> = 0;

  virtual auto getHDivNodesInternal(indexType meshID, indexType order)
      -> std::vector<GenericNodes *> = 0;

  virtual void setHDivShapes(indexType meshid,
                             indexType order, NodeTypes type) = 0;

  virtual void getHDivDofs(std::vector<DegreeOfFreedom *> &Dofs,
                           indexType meshID, indexType order,
                           NodeTypes type) = 0;

  virtual void getHDivShapes(indexType order,
                             Types::Matrix2X<prec> &shape,
                             Types::VectorX<prec> &dshape, prec xi,
                             prec eta) = 0;

  virtual auto getHDivShapes(indexType order,
                             IntegrationPoint &IntegrationPt) -> HDivShapes = 0;
};

} // namespace HierAMuS::Geometry
