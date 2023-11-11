// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "geometry/GeometryBaseData.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

#include <geometry/GeometryData.h>

#include "geometry/Volumes/LinearBrickRuntime.h"
#include <geometry/Edges/EdgesData.h>
#include <geometry/Faces/FacesData.h>
#include <geometry/VertexData.h>
#include <geometry/Volumes/LinearBrickData.h>

#include <types/MatrixTypes.h>

#include <vector>

#include <iomanip>

#include "plot/vtkplotClass.h"

#include "shapefunctions/LobattoShapes.h"
#include "shapefunctions/IntegrationsPoints/helperClasses/IntegrationPointsManagement.h"

#include <vtkCellType.h>
namespace HierAMuS::Geometry {

LinearBrickData::LinearBrickData() {}

LinearBrickData::~LinearBrickData() = default;

auto LinearBrickData::getType() -> const GeometryTypes & {
  return HierAMuS::Geometry::LinearBrickData::type;
}

auto LinearBrickData::getRuntimeObject(GeometryData &geoData)
    -> std::shared_ptr<VolumesRuntime> {
  return std::make_shared<LinearBrickRuntime>(geoData,
                                              *this);
}

void LinearBrickData::setAllNodeBoundaryConditionMeshId(indexType meshId,
                                                        indexType dof) {
  GeometryBaseData::setAllNodeBoundaryConditionMeshId(meshId, dof);

  for (auto &i : m_faces_pointers) {
    i->setAllNodeBoundaryConditionMeshId(meshId, dof);
  }
}

auto LinearBrickData::getIntegrationPoints(indexType elementId)
    -> IntegrationPoints {
  auto points = IntegrationPointsManagement::getIntegrationsPoints(elementId);
  points.setType(IntegrationType::Gauss3D);
  return points;
}

void LinearBrickData::setH1Shapes(indexType meshId, indexType order,
                                  NodeTypes type) {
  for (auto i = 0; i < m_numberOfFaces; i++) {
    m_faces_pointers[i]->setH1Shapes(meshId, order, type);
  }
  if (order > 1) {
    this->setH1ShapesInternal(meshId, order, type);
  }
}

void LinearBrickData::setH1ShapesInternal(indexType meshId, indexType order,
                                          NodeTypes type) {
  if (order > 1) {
    indexType num = order - 1;
    num *= num * num;
    this->setNodeSet(meshId, num, type);
  }
}

void LinearBrickData::checkUpdateElement(EquationHandler &eqHandler,
                                         GeometryData &geoData) {

  static const std::vector<indexType> lv1({0, 0, 1, 2, 3, 4});
  static const std::vector<indexType> lv2({1, 1, 2, 3, 0, 5});
  static const std::vector<indexType> lv3({2, 5, 6, 7, 4, 6});
  static const std::vector<indexType> lv4({3, 4, 5, 6, 7, 7});

  for (auto i = 0; i < 6; ++i) {
    if (this->m_faces[i] == -1) {
      auto &V1 = geoData.getVertexData(this->m_verts[lv1[i]]);
      auto &V2 = geoData.getVertexData(this->m_verts[lv2[i]]);
      auto &V3 = geoData.getVertexData(this->m_verts[lv3[i]]);
      auto &V4 = geoData.getVertexData(this->m_verts[lv4[i]]);

      auto faceNums = V1.getConnectedFaces();
      bool faceExists = false;
      bool search = true;
      if (faceNums.empty()) {
        search = false;
      }
      indexType pos = 0;
      while (search) {
        auto face = geoData.getFaceData(faceNums[pos]);
        if (face->hasVertices(V1.getId(), V2.getId(), V3.getId())) {
          faceExists = true;
          search = false;
          this->m_faces[i] = face->getId();
        } else {
          ++pos;
          if (pos >= static_cast<indexType>(faceNums.size())) {
            search = false;
          }
        }
      }
      if (!faceExists) {
        auto facenum = geoData.requestNewGeometryObject(
            eqHandler, GeometryTypes::LinearQuadrilateral);
        auto face = geoData.getFaceData(facenum);
        std::vector<indexType> verts(
            {V1.getId(), V2.getId(), V3.getId(), V4.getId()});
        face->setVerts(geoData, verts);
        this->m_faces[i] = face->getId();
        face->checkUpdateElement(eqHandler, geoData);
      }
    }
  }

  // Edges
  static const std::vector<indexType> le1({0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7});
  static const std::vector<indexType> le2({1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 7, 4});
  for (auto i = 0; i < 12; ++i) {
    if (this->m_edges[i] == -1) {
      auto &V1 = geoData.getVertexData(this->m_verts[le1[i]]);
      auto &V2 = geoData.getVertexData(this->m_verts[le2[i]]);

      bool search = true;
      bool edgeExists = false;
      auto edgeNums = V1.getConnectedEdges();
      if (edgeNums.empty()) {
        search = false;
      }
      indexType pos = 0;
      while (search) {
        auto &edge = geoData.getEdgeData(edgeNums[pos]);
        if (edge.hasVertices(V1.getId(), V2.getId())) {
          edgeExists = true;
          search = false;
          this->m_edges[i] = edge.getId();
        } else {
          pos++;
          if (pos >= static_cast<indexType>(edgeNums.size())) {
            search = false;
          }
        }
      }
      if (!edgeExists) {
        std::cout << "something went wrong" << std::endl;
        auto edgeNum = geoData.requestNewGeometryObject(
            eqHandler, GeometryTypes::LinearEdge);
        auto &edge = geoData.getEdgeData(edgeNum);
        std::vector<indexType> verts = {V1.getId(), V2.getId()};
        edge.setVerts(geoData, verts);
        this->m_edges[i] = edgeNum;
      }
    }
  }
}

const GeometryTypes LinearBrickData::type = GeometryTypes::LinearBrick;
} // namespace HierAMuS::Geometry
