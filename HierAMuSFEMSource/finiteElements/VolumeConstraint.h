// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "Volume.h"
#include "geometry/Faces.h"
#include "geometry/Volumes.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

namespace HierAMuS::FiniteElement {
class VolumeConstraint : public Volume {
  using ptrCol = PointerCollection;

public:
  VolumeConstraint() = default;
  ~VolumeConstraint() override;

  void setVerts(std::vector<indexType> &vertsIn) override;

  auto getVertexCoordinates(ptrCol &pointers) -> Types::Vector3<prec>;
  auto getVolumeCoordinates(ptrCol &pointers, IntegrationPoint &IntegrationPt)
      -> Types::Vector3<prec>;
  
  void setVertexNodes(PointerCollection &pointers, indexType meshId);
  void getVertexDofs(PointerCollection &pointers,
                     std::vector<DegreeOfFreedom *> &Dofs, indexType meshID);

  auto getVertexNodes(PointerCollection &pointers, indexType meshId)
      -> std::vector<GenericNodes *>;
  
private:
  indexType sharedVertex;
};
} // namespace HierAMuS
