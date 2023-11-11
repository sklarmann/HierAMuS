// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "geometry/Edges/EdgesRuntime.h"
#include "geometry/Faces/FacesRuntime.h"

#include "geometry/GeometryShape.h"
#include "geometry/VertexRuntime.h"
#include <geometry/GeometryTypes.h>
#include <types/MatrixTypes.h>

#include <vector>

namespace HierAMuS::Geometry {
class GeometryData;
class VolumesData;
class FacesRuntime;

// General Interface to switch between fixed and variable size
template <bool, indexType numVertices, class TData,
          class self>
class EdgesRuntimeDataInterfaceVariables : public EdgesRuntime {};

// Specialization for fixed array size
template <indexType numVertices, class TData, class self>
class EdgesRuntimeDataInterfaceVariables<true, numVertices, TData, self>
    : public EdgesRuntime {
public:
  EdgesRuntimeDataInterfaceVariables(GeometryData &geoData,
                                     TData &data_element)
      : EdgesRuntime(data_element), m_Edge_Data_Element(data_element) {
    for (auto i = 0; i < m_numberOfVerts; ++i) {
      m_Vertices[i] =
          geoData.getVertexRuntimePointer(data_element.getVertexNumber(i));
    }
  }

protected:
  static constexpr indexType m_numberOfVerts = numVertices;

  TData &m_Edge_Data_Element;
  std::array<VertexRuntime *, m_numberOfVerts> m_Vertices;
};

// Specialization for dynamic array size
template <indexType numVertices, class TData, class self>
class EdgesRuntimeDataInterfaceVariables<false, numVertices, TData,
                                         self> : public EdgesRuntime {
  using baseBool = typename std::conditional_t<
          numVertices <= 0,
          std::integral_constant<bool, false>,
          std::integral_constant<bool, true>>;
public:
  EdgesRuntimeDataInterfaceVariables(GeometryData &geoData,
                                     TData &data_element)
      : EdgesRuntime(data_element), m_Edges_Data_Element(data_element) {
  
    m_numberOfVerts = data_element.getNumberOfVerts();

    m_Vertices.resize(m_numberOfVerts);
    for (auto i = 0; i < m_numberOfVerts; ++i) {
      m_Vertices[i] =
          geoData.getVertexRuntimePointer(data_element.getVertexNumber(i));
    }
  }

protected:
  indexType m_numberOfVerts;

  TData &m_Edges_Data_Element;
  std::vector<VertexRuntime *> m_Vertices;
};

template <indexType numVertices, class TData, class self>
class EdgesRuntimeDataInterface
    : public std::conditional_t<numVertices <= 0,
                                EdgesRuntimeDataInterfaceVariables<
                                    false, numVertices, TData, self>,
                                EdgesRuntimeDataInterfaceVariables<
                                    true, numVertices, TData, self>> {
  using baseBool = typename std::conditional_t<
          numVertices <= 0,
          std::integral_constant<bool, false>,
          std::integral_constant<bool, true>>;
public:
  EdgesRuntimeDataInterface(GeometryData &geoData, TData &data_element)
      : EdgesRuntimeDataInterfaceVariables<baseBool::value, numVertices, TData, self>(geoData,data_element) {

  };
  ~EdgesRuntimeDataInterface() = default;

  /**
   * @brief Get a vector with pointers of type Geometry::VertexRume to the
   * vertices of the geometric element.
   *
   * @param edgesOut[out], a vector of size 12 with pointers to the edges of the
   * geometric element.
   */
  void getVertices(std::vector<VertexRuntime *> &verticesOut) override {
    verticesOut =
        std::vector<VertexRuntime *>(this->m_Vertices.begin(), this->m_Vertices.end());
  };
  /**
   * @brief Get a pointer to the local vertex local_number of the quadrilateral
   * face.
   *
   * @param pointers object containing the pointers to global data.
   * @param local_number local number of the vertex.
   * @return Geometry::Vertex* pointer to the vertex.
   */
  auto getVertex(indexType local_number) -> Geometry::VertexRuntime & override {
    return *this->m_Vertices[local_number];
  };

  /**
   * @brief Returns the global number of the local vertex with number
   * localNumber.
   *
   * @param[in] localNumber, the local number of the vertex.
   * @return indexType, the global number of the vertex.
   */
  auto getVertexNumber(indexType localNumber) -> indexType override {
    return this->m_Vertices[localNumber]->getId();
  };
  /**
   * @brief Get the global vertex numbers of the geometric element.
   *
   * @param vertsOut[out], a vector with the global numbers of the vertices of
   * type indexType and of size 8.
   */
  auto getVertexNumbers() -> std::vector<indexType> override {
    std::vector<indexType> nums(this->m_numberOfVerts);
    nums.resize(this->m_numberOfVerts);
    for (auto i = 0; i < this->m_numberOfVerts; ++i) {
      nums[i] = this->m_Vertices[i]->getId();
    }
    return nums;
  };

  /**
   * @brief Get the Number Of Vertices of the geometric object, return 8
   *
   * @return indexType, returns 8.
   */
  auto getNumberOfVerts() -> indexType override { return this->m_numberOfVerts; };

  /**
   * @brief Get coordinate at integration point.
   *
   * @param pointers, pointer to global data.
   * @param IntPoint, integration point to evaluate coordinates at.
   *
   * @return Vector3 with the x,y,z coordinates.
   */
  auto getCoordinates(IntegrationPoint &IntPoint)
      -> Types::Vector3<prec> override {
    auto shp = static_cast<self *>(this)->getH1Shapes(1, IntPoint);

    Types::Vector3<prec> coor(Types::Vector3<prec>::Zero());
    for (indexType i = 0; i < this->m_numberOfVerts; ++i) {
      coor += shp.shapes(i) * this->m_Vertices[i]->getCoordinates();
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
      -> prec override {
    
    Types::Vector3<prec> temp =
        this->m_Vertices[this->m_numberOfVerts - 1]->getCoordinates() - this->m_Vertices[0]->getCoordinates();
    prec jacobi = temp.stableNorm();
    return jacobi;
  };
  /**
   * @brief Changes orientation of the face (reverts the face normal vector).
   */
  void flip() override {
    std::reverse(this->m_Vertices.begin() + 1, this->m_Vertices.end());
  }


  /** @brief Check if the face has the three vertices.
   *
   * @param[in] v1 Global number of a corner vertex.
   * @param[in] v2 Global number of a corner vertex.
   * @param[in] v3 Global number of a corner vertex.
   * @return true If the face has the three vertices.
   * @return false If the face does not have the three vertices.
   */
  auto hasVertices(indexType v1, indexType v2) -> bool override {
    bool t1 = false, t2 = false;
    for (auto i = 0; i < this->m_numberOfVerts; ++i) {
      if (this->m_Vertices[i]->getId() == v1) {
        t1 = true;
      }
      if (this->m_Vertices[i]->getId() == v2) {
        t2 = true;
      }
    }

    return (t1 && t2);
  };

  void set_geometry_pointers(GeometryData &geoData) override {
    for (auto i = 0; i < this->m_numberOfVerts; ++i) {
      this->m_Vertices[i] = geoData.getVertexRuntimePointer(
          this->m_Edge_Data_Element.getVertexNumber(i));
    }
  };

protected:
};

} // namespace HierAMuS::Geometry
