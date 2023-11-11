// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "VolumeConstraint.h"
#include "MatrixTypes.h"
#include "geometry/GeometryBaseData.h"
#include "geometry/GeometryData.h"
#include "GenericNodes.h"
#include <vector>
#include "geometry/Volumes/VolumesRuntime.h"
#include "geometry/Volumes/VolumesH1Interface.h"

namespace HierAMuS::FiniteElement {
VolumeConstraint::~VolumeConstraint() = default;

void VolumeConstraint::set_pointers(PointerCollection &pointers) {
  m_volume_pointer = pointers.getGeometryData()->getVolumeData(m_volume);
  m_vertex_pointer = &pointers.getGeometryData()->getVertexData(m_sharedVertex);

  m_Vertex =
      pointers.getGeometryData()->getVertexRuntimePointer(m_sharedVertex);
  m_Volume = pointers.getGeometryData()->getVolumeRuntime(m_volume);
}

void VolumeConstraint::setVerts(std::vector<indexType> &vertsIn)
{
  m_sharedVertex = vertsIn[0];
}

auto VolumeConstraint::getVertexCoordinates(ptrCol &pointers)
    -> Types::Vector3<prec> {
  return pointers.getGeometryData()->getVertexData(m_sharedVertex).getCoordinates();
}

auto VolumeConstraint::getVolumeCoordinates(ptrCol &pointers,
                                            IntegrationPoint &IntegrationPt)
    -> Types::Vector3<prec> {
  return m_Volume->getCoordinates(IntegrationPt);
}

void VolumeConstraint::setVertexNodes(PointerCollection &pointers,
                                      indexType meshId)
{
  auto &Vert = pointers.getGeometryData()->getVertexData(m_sharedVertex);
  Vert.setNodeSet(meshId, 1, NodeTypes::displacement);
}

void VolumeConstraint::getVertexDofs(PointerCollection &pointers,
                                     std::vector<DegreeOfFreedom *> &Dofs,
                                     indexType meshID) {
  auto &vertex = pointers.getGeometryData()->getVertexData(this->m_sharedVertex);
  std::vector<GenericNodes *> nodeVector;
  vertex.getNodes(nodeVector, meshID);
  for (auto &node : nodeVector) {
    auto addDofs = node->getDegreesOfFreedom();
    Dofs.insert(Dofs.end(), addDofs.begin(), addDofs.end());
  }
}

auto VolumeConstraint::getVertexNodes(PointerCollection &pointers,
                                      indexType meshId)
    -> std::vector<GenericNodes *> {

  auto &V = pointers.getGeometryData()->getVertexData(m_sharedVertex);
  return V.getNodesOfSet(meshId);
}

void VolumeConstraint::setVolume(indexType volumeIn) {
  this->m_volume = volumeIn;
}

auto VolumeConstraint::getJacobian(ptrCol &pointers,
                                   IntegrationPoint &IntegrationPt)
    -> Types::Matrix33<prec> {

  return m_Volume->getJacobian(IntegrationPt);
}

auto VolumeConstraint::getH1Shapes(ptrCol &pointers, indexType order,
                                   Types::Matrix33<prec> &jacobi,
                                   IntegrationPoint &IntegrationPt)
    -> Geometry::H1Shapes {
  auto shapes = m_Volume->getH1Volume()->getH1Shapes(order, IntegrationPt);
  Types::Matrix33<prec> jacobiInv;
  jacobiInv = jacobi.inverse().transpose();

  shapes.shapeDeriv = jacobiInv * shapes.shapeDeriv;
  return shapes;
}

auto VolumeConstraint::getIntegrationPoints(ptrCol &pointers)
    -> IntegrationPoints {
  auto temp = pointers.getGeometryData()->getVolumeData(this->m_volume);
  return temp->getIntegrationPoints(this->m_id);
  ;
}

void VolumeConstraint::setH1Shapes(ptrCol &pointers, indexType meshid,
                                   indexType order) {
  auto tempVol = pointers.getGeometryData()->getVolumeData(this->m_volume);
  tempVol->setH1Shapes(meshid, order, NodeTypes::displacement);
}

void VolumeConstraint::getH1Dofs(ptrCol &pointers,
                                 std::vector<DegreeOfFreedom *> &Dofs,
                                 indexType meshID, indexType order) {
  m_Volume->getH1Volume()->getH1Dofs(Dofs, meshID, order);
}

auto VolumeConstraint::getH1Nodes(ptrCol &pointers, indexType meshID,
                                  indexType order)
    -> std::vector<GenericNodes *> {
  return m_Volume->getH1Volume()->getH1Nodes(meshID, order);
}

} // namespace HierAMuS
