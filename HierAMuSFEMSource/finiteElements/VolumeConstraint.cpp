// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "VolumeConstraint.h"
#include "MatrixTypes.h"
#include "geometry/Base.h"
#include <vector>

namespace HierAMuS::FiniteElement {
VolumeConstraint::~VolumeConstraint() = default;

void VolumeConstraint::setVerts(std::vector<indexType> &vertsIn)
{
  sharedVertex = vertsIn[0];
}

auto VolumeConstraint::getVertexCoordinates(ptrCol &pointers)
    -> Types::Vector3<prec> {
  return pointers.getGeometryData()->getVertex(sharedVertex).getCoordinates();
}

auto VolumeConstraint::getVolumeCoordinates(ptrCol &pointers,
                                            IntegrationPoint &IntegrationPt)
    -> Types::Vector3<prec> {
  return pointers.getGeometryData()->getVolume(volume)->getCoordinates(pointers,IntegrationPt);
}

void VolumeConstraint::setVertexNodes(PointerCollection &pointers,
                                      indexType meshId)
{
  auto &Vert = pointers.getGeometryData()->getVertex(sharedVertex);
  Vert.setNodeSet(pointers, meshId, 1, NodeTypes::displacement);
}

void VolumeConstraint::getVertexDofs(PointerCollection &pointers,
                                     std::vector<DegreeOfFreedom *> &Dofs,
                                     indexType meshID) {
  auto &vertex = pointers.getGeometryData()->getVertex(this->sharedVertex);
  std::vector<GenericNodes *> nodeVector;
  vertex.getNodes(pointers, nodeVector, meshID);
  for (auto &node : nodeVector) {
    auto addDofs = node->getDegreesOfFreedom(pointers);
    Dofs.insert(Dofs.end(), addDofs.begin(), addDofs.end());
  }
}

auto VolumeConstraint::getVertexNodes(PointerCollection &pointers,
                                      indexType meshId)
    -> std::vector<GenericNodes *> {

  auto &V = pointers.getGeometryData()->getVertex(sharedVertex);
  return V.getNodesOfSet(pointers, meshId);
}

} // namespace HierAMuS
