// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "datatypes.h"
#include "Nodetypes.h"
#include "MeshIdNodeList.h"
#include "geometry/GeometryShape.h"
#include "geometry/Edges/EdgesRuntime.h"


namespace HierAMuS {
class DegreeOfFreedom;
class GenericNodes;
class IntegrationPoint;
}

namespace HierAMuS::Geometry {
class EdgesRuntime;

class EdgeH1ShapesInterface {
private:

public:
  //EdgeH1ShapesInterface() = default;
  ~EdgeH1ShapesInterface() = default;

  virtual auto getH1Edge() -> EdgeH1ShapesInterface * = 0;

  virtual void setH1Shapes(indexType meshId, indexType order,
                           NodeTypes type) = 0;
  virtual void setH1ShapesInternal(indexType meshId, indexType order,
                                   NodeTypes type) = 0;

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

  virtual auto getH1Shapes(indexType order, IntegrationPoint &integration_point)
      -> H1Shapes = 0;

  virtual auto getH1ShapesInternal(indexType order,
                                   IntegrationPoint &integration_point)
      -> H1Shapes = 0;

};
} // namespace HierAMuS::Geometry
