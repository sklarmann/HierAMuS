// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <geometry/Volumes/VolumesRuntime.h>

#include "geometry/VertexRuntime.h"
#include "geometry/GeometryShape.h"
#include <geometry/GeometryTypes.h>
#include <types/MatrixTypes.h>

#include <vector>

namespace HierAMuS::Geometry {
class GeometryData;
class VolumesData;
class FacesRuntime;

template <indexType numVertices, indexType numEdges, indexType numFaces, class TData, class self>
class VolumesRuntimeDataInterface : public VolumesRuntime {
public:
  VolumesRuntimeDataInterface(GeometryData &geoData, TData &data_element)
      : VolumesRuntime(data_element), m_Volume_Data_Element(data_element) {

    for (auto i = 0; i < m_numberOfVerts; ++i) {
      m_Vertices[i] = geoData.getVertexRuntimePointer(
          data_element.getVertexNumber(i));
    }
    for (auto i = 0; i < m_numberOfEdges; ++i) {
      m_Edges[i] =
          geoData.getEdgeRuntime(data_element.getEdgeNumber(i));
    }
    for (auto i = 0; i < m_numberOfFaces; ++i) {
      m_Faces[i] =
          geoData.getFaceRuntime(data_element.getFaceNumber(i));
    }
  };
  ~VolumesRuntimeDataInterface() = default;

  /**
   * @brief Get a vector with pointers of type Geometry::VertexRume to the vertices of
   * the geometric element.
   *
   * @param edgesOut[out], a vector of size 12 with pointers to the edges of the
   * geometric element.
   */
  void getVertices(std::vector<VertexRuntime *> &verticesOut) override {
    verticesOut.resize(m_numberOfVerts);
    for (auto i = 0; i < m_numberOfVerts; ++i) {
      verticesOut[i] = m_Vertices[i];
    }
  };
  /**
   * @brief Get a vector with pointers of type Geometry::Base to the edges of
   * the geometric element.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param edgesOut[out], a vector of size 12 with pointers to the edges of the
   * geometric element.
   */
  void getEdges(std::vector<EdgesRuntime *> &edgesOut) override{
    edgesOut.resize(m_numberOfEdges);
    for (auto i = 0; i < m_numberOfEdges; ++i) {
      edgesOut[i] = m_Edges[i].get();
    }
  };
  /**
   * @brief Get a vector with pointers of type Geometry::Base to the faces of
   * the geometric element.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param facesOut[out], a vector of size 6 with pointers to the faces of the
   * geometric element.
   */
  void getFaces(std::vector<FacesRuntime *> &facesOut) override {
    facesOut.resize(m_numberOfFaces);
    for (auto i = 0; i < m_numberOfFaces; ++i) {
      facesOut[i] = m_Faces[i].get();
    }
  };
  /**
   * @brief Get the global vertex numbers of the geometric element.
   *
   * @param vertsOut[out], a vector with the global numbers of the vertices of
   * type indexType and of size 8.
   */
  auto getVertexNumbers() -> std::vector<indexType> override {
    std::vector<indexType> nums(m_numberOfVerts);
    nums.resize(m_numberOfVerts);
    for (auto i = 0; i < m_numberOfVerts; ++i) {
      nums[i] = m_Vertices[i]->getId();
    }
    return nums;
  };

  /**
   * @brief Get the Number Of Vertices of the geometric object, return 8
   *
   * @return indexType, returns 8.
   */
  auto getNumberOfVerts() -> indexType override { return m_numberOfVerts; };
  /**
   * @brief Get the Number of Edges of the element, returns 12.
   *
   * @return indexType, returns 12.
   */
  auto getNumberOfEdges() -> indexType override { return m_numberOfEdges; };
  /**
   * @brief Get the Number Of Faces object
   *
   * @return indexType Number of faces as an indexType
   */
  auto getNumberOfFaces() -> indexType override { return m_numberOfFaces; };
  /**
   * @brief Get coordinate at integration point.
   *
   * @param pointers, pointer to global data.
   * @param IntPoint, integration point to evaluate coordinates at.
   *
   * @return Vector3 with the x,y,z coordinates.
   */
  auto getCoordinates(IntegrationPoint &IntPoint)
      -> Types::Vector3<prec> override{
    auto shp = static_cast<self *>(this)->getH1Shapes(
        1, IntPoint);

    Types::Vector3<prec> coor(Types::Vector3<prec>::Zero());
    for (indexType i = 0; i < m_numberOfVerts; ++i) {
      coor += shp.shapes(i) * m_Vertices[i]->getCoordinates();
    }

    return coor;
  };
  /**
   * @brief Get the Jacobian matrix at the current IntegrationPoint. New version
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param IntegrationPt[in], current IntegrationPoint.
   * @return Types::MatrixXX<prec>, Jacobian matrix at the current
   * IntegrationPoint.
   */
  auto getJacobian(IntegrationPoint &IntegrationPt)
      -> Types::Matrix33<prec> override {
    Types::Matrix33<prec> jacobi(Types::Matrix33<prec>::Zero());
    auto shp = static_cast<self *>(this)->getH1Shapes(
        1, IntegrationPt);
    for (auto i = 0; i < m_numberOfVerts; ++i) {
      for (auto j = 0; j < 3; j++) {
        jacobi.block(0, j, 3, 1) +=
            shp.shapeDeriv(j, i) * m_Vertices[i]->getCoordinates();
      }
    }
    return jacobi;
  };
  

protected:
  static constexpr indexType m_numberOfVerts = numVertices;
  static constexpr indexType m_numberOfEdges = numEdges;
  static constexpr indexType m_numberOfFaces = numFaces;

  TData &m_Volume_Data_Element;
  std::array<VertexRuntime *, m_numberOfVerts> m_Vertices;
  std::array<std::shared_ptr<EdgesRuntime>, m_numberOfEdges> m_Edges;
  std::array<std::shared_ptr<FacesRuntime>, m_numberOfFaces> m_Faces;
};

} // namespace HierAMuS::Geometry
