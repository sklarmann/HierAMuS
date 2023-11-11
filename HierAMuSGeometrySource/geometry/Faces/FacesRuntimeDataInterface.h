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
template <bool, indexType numVertices, indexType numEdges, class TData,
          class self>
class FacesRuntimeDataInterfaceVariables : public FacesRuntime {};

// Specialization for fixed array size
template <indexType numVertices, indexType numEdges, class TData, class self>
class FacesRuntimeDataInterfaceVariables<true, numVertices, numEdges, TData,
                                         self> : public FacesRuntime {
public:
  FacesRuntimeDataInterfaceVariables(GeometryData &geoData,
                                     TData &data_element)
      : FacesRuntime(data_element), m_Face_Data_Element(data_element) {
    for (auto i = 0; i < m_numberOfVerts; ++i) {
      m_Vertices[i] =
          geoData.getVertexRuntimePointer(data_element.getVertexNumber(i));
    }
    for (auto i = 0; i < m_numberOfEdges; ++i) {
      m_Edges[i] = geoData.getEdgeRuntime(data_element.getEdgeNumber(i));
    }
  }

protected:
  static constexpr indexType m_numberOfVerts = numVertices;
  static constexpr indexType m_numberOfEdges = numEdges;

  TData &m_Face_Data_Element;
  std::array<VertexRuntime *, m_numberOfVerts> m_Vertices;
  std::array<std::shared_ptr<EdgesRuntime>, m_numberOfVerts> m_Edges;
};

// Specialization for dynamic array size
template <indexType numVertices, indexType numEdges, class TData, class self>
class FacesRuntimeDataInterfaceVariables<false, numVertices, numEdges, TData,
                                         self> : public FacesRuntime {
public:
  FacesRuntimeDataInterfaceVariables(GeometryData &geoData,
                                     TData &data_element)
      : FacesRuntime(data_element),
        m_Face_Data_Element(data_element) {
  
    m_numberOfVerts = data_element.getNumberOfVerts();
    m_numberOfEdges = data_element.getNumberOfEdges();

    m_Vertices.resize(m_numberOfVerts);
    for (auto i = 0; i < m_numberOfVerts; ++i) {
      m_Vertices[i] =
          geoData.getVertexRuntimePointer(data_element.getVertexNumber(i));
    }
    m_Edges.resize(m_numberOfEdges);
    for (auto i = 0; i < m_numberOfEdges; ++i) {
      m_Edges[i] = geoData.getEdgeRuntime(data_element.getEdgeNumber(i));
    }

  }

protected:
  indexType m_numberOfVerts;
  indexType m_numberOfEdges;

  TData &m_Face_Data_Element;
  std::vector<VertexRuntime *> m_Vertices;
  std::vector<std::shared_ptr<EdgesRuntime>> m_Edges;
};

template <indexType numVertices, indexType numEdges, class TData, class self>
class FacesRuntimeDataInterface
    : public std::conditional_t<numVertices <= 0 || numEdges <= 0,
                                FacesRuntimeDataInterfaceVariables<
                                    false, numVertices, numEdges, TData, self>,
                                FacesRuntimeDataInterfaceVariables<
                                    true, numVertices, numEdges, TData, self>> {
  using baseBool = typename std::conditional_t<
          numVertices <= 0 || numEdges <= 0,
          std::integral_constant<bool, false>,
          std::integral_constant<bool, true>>;
public:
  FacesRuntimeDataInterface(GeometryData &geoData, TData &data_element)
      : FacesRuntimeDataInterfaceVariables<baseBool::value,numVertices,numEdges, TData, self>(geoData,data_element) {

  };
  ~FacesRuntimeDataInterface() = default;

  /**
   * @brief Get a vector with pointers of type Geometry::VertexRume to the
   * vertices of the geometric element.
   *
   * @param edgesOut[out], a vector of size 12 with pointers to the edges of the
   * geometric element.
   */
  void getVertices(std::vector<VertexRuntime *> &verticesOut) override {
    verticesOut.resize(this->m_numberOfVerts);
    for (auto i = 0; i < this->m_numberOfVerts; ++i) {
      verticesOut[i] = this->m_Vertices[i];
    }
  };
  /**
   * @brief Get a pointer to the local vertex local_number of the quadrilateral
   * face.
   *
   * @param pointers object containtig the pointers to global data.
   * @param local_number local number of the vertex.
   * @return Geometry::Vertex* pointer to the vertex.
   */
  auto getVertex(indexType local_number) -> Geometry::VertexRuntime * override {
    return this->m_Vertices[local_number];
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
   * @brief Get a vector with pointers of type Geometry::Base to the edges of
   * the geometric element.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param edgesOut[out], a vector of size 12 with pointers to the edges of the
   * geometric element.
   */
  void getEdges(std::vector<EdgesRuntime *> &edgesOut) override {
    edgesOut.resize(this->m_numberOfEdges);
    for (auto i = 0; i < this->m_numberOfEdges; ++i) {
      edgesOut[i] = this->m_Edges[i].get();
    }
  };
  /**
   * @brief Get a pointer to the local edge local_number of the quadrilateral
   * face.
   *
   * @param pointers object containtig the pointers to global data.
   * @param local_number local number of the edge.
   * @return Geometry::Vertex* pointer to the edge.
   */
  auto getEdge(indexType local_number) -> Geometry::EdgesRuntime * override {
    return this->m_Edges[local_number].get();
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
   * @brief Get the global edge numbers of the object.
   *
   * @param[out] edgesOut, vector of 4 integers of type indexType with global
   * edge numbers of the element.
   */
  auto getEdgeNumbers() -> std::vector<indexType> override {
    std::vector<indexType> edgeNumber(this->m_numberOfEdges);
    for (auto i = 0; i < this->m_numberOfEdges; ++i) {
      edgeNumber[i] = this->m_Edges[i]->getId();
    }
    return edgeNumber;
  };
  /**
   * @brief Get the global edge numbers of the object.
   *
   * @param[out] edgesOut, vector of 4 integers of type indexType with global
   * edge numbers of the element.
   */
  auto getEdgeNumber(indexType local_number) -> indexType override {
    return this->m_Edges[local_number]->getId();
  };

  /**
   * @brief Get the Number Of Vertices of the geometric object, return 8
   *
   * @return indexType, returns 8.
   */
  auto getNumberOfVerts() -> indexType override { return this->m_numberOfVerts; };
  /**
   * @brief Get the Number of Edges of the element, returns 12.
   *
   * @return indexType, returns 12.
   */
  auto getNumberOfEdges() -> indexType override { return this->m_numberOfEdges; };

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
      -> Types::Matrix22<prec> override {
    Types::Matrix22<prec> jacobi(Types::Matrix22<prec>::Zero());
    auto shp = static_cast<self *>(this)->getH1Shapes(1, IntegrationPt);
    for (auto i = 0; i < this->m_numberOfVerts; ++i) {
      for (auto j = 0; j < 2; j++) {
        jacobi.block(0, j, 2, 1) +=
            shp.shapeDeriv(j, i) *
            this->m_Vertices[i]->getCoordinates().block(0, 0, 2, 1);
      }
    }
    return jacobi;
  };
  /**
   * @brief Changes orientation of the face (reverts the face normal vector).
   */
  void flip() override {
    std::reverse(this->m_Vertices.begin() + 1, this->m_Vertices.end());
    std::reverse(this->m_Edges.begin(), this->m_Edges.end());
  }
  /**
   * @brief Rotates the face clockwise n times.
   *
   * @param [in] n Number of rotations.
   */
  void rotate(indexType n) override {
    while (n >= this->m_numberOfVerts) {
      n -= this->m_numberOfVerts;
    }
    while (n < this->m_numberOfVerts) {
      n += this->m_numberOfVerts;
    }
    std::rotate(this->m_Vertices.begin(), this->m_Vertices.begin() + n,
                this->m_Vertices.end());
    std::rotate(this->m_Edges.begin(), this->m_Edges.begin() + n,
                this->m_Edges.end());
  }

  /**
   * @brief Computes the mean coordinate of the face element.
   *
   * @param [in] pointers Pointer collection to global data.
   *
   * @return Mean coordinate of the face element.
   */
  auto computeMeanCoordinate()
      -> Types::Vector3<prec> override {

    Types::Vector3<prec> meanCoor = Types::Vector3<prec>::Zero();

    for (auto &i : this->m_Vertices) {
      meanCoor += i->getCoordinates();
    }
    meanCoor /= prec(this->m_numberOfVerts);
    return meanCoor;
  }

  /** @brief Check if the face has the three vertices.
   *
   * @param[in] v1 Global number of a corner vertex.
   * @param[in] v2 Global number of a corner vertex.
   * @param[in] v3 Global number of a corner vertex.
   * @return true If the face has the three vertices.
   * @return false If the face does not have the three vertices.
   */
  auto hasVertices(indexType v1, indexType v2, indexType v3) -> bool override {
    bool t1 = false, t2 = false, t3 = false;
    for (auto i = 0; i < this->m_numberOfVerts; ++i) {
      if (this->m_Vertices[i]->getId() == v1) {
        t1 = true;
      }
      if (this->m_Vertices[i]->getId() == v2) {
        t2 = true;
      }
      if (this->m_Vertices[i]->getId() == v3) {
        t3 = true;
      }
    }

    return (t1 && t2 && t3);
  };

protected:
};

} // namespace HierAMuS::Geometry
