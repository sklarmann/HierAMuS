// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "geometry/Faces/FacesData.h"


#include <type_traits>

#include "HelperFunctions.h"
#include "geometry/GeometryShape.h"
#include "geometry/VertexRuntime.h"
#include <geometry/GeometryTypes.h>
#include <types/MatrixTypes.h>

#include <algorithm>
#include <iterator>
#include <vector>

namespace HierAMuS::Geometry {
class GeometryData;
class FacesRuntime;

template <bool, indexType numVertices>
class EdgesDataInterfaceVariables : public EdgesData {};

// Static sized storage
template <indexType numVertices>
class EdgesDataInterfaceVariables<true, numVertices> : public EdgesData {
public:
  EdgesDataInterfaceVariables()
      : EdgesData(), m_verts(initArray<numVertices, indexType, -1>()){};

protected:
  inline static const indexType m_numberOfVerts = numVertices;

  std::array<indexType, numVertices> m_verts;

  std::array<VertexData *, numVertices> m_verts_pointers;
};

// Dynamic sized storage
template <indexType numVertices>
class EdgesDataInterfaceVariables<false, numVertices>
    : public EdgesData {
public:
  EdgesDataInterfaceVariables()
      : EdgesData(), m_numberOfVerts(0) {};

  void setNumberOfVertices(indexType numberOfVertices) {
    m_numberOfVerts = numberOfVertices;
    m_verts.resize(numberOfVertices);
    m_verts_pointers.resize(numberOfVertices);
    std::fill(m_verts.begin(), m_verts.end(), -1);
  }

protected:
  indexType m_numberOfVerts;

  std::vector<indexType> m_verts;

  std::vector<VertexData *> m_verts_pointers;
};

template <indexType numVertices, class self>
class EdgesDataInterface
    : public std::conditional_t<
          numVertices <= 0,
          EdgesDataInterfaceVariables<false, numVertices>,
          EdgesDataInterfaceVariables<true, numVertices>>
{
  using baseBool = typename std::conditional_t<
          numVertices <= 0,
          std::integral_constant<bool, false>,
          std::integral_constant<bool, true>>;
public:
  EdgesDataInterface() : EdgesDataInterfaceVariables<baseBool::value,numVertices>(){};
  ~EdgesDataInterface() = default;
  /**
   * @brief Set the Vertex numbers of the geometric element.
   * They need to be given according to the numbering in the picture above.
   *
   * @param vertsIn[in], a vector with the global numbers of the vertices of
   * type indexType and of size 8.
   */
  void setVerts(GeometryData &geoData,
                std::vector<indexType> &vertsIn) override {
    std::copy_n(vertsIn.begin(), this->m_numberOfVerts, this->m_verts.begin());
    for (auto &i : this->m_verts) {
      auto &V = geoData.getVertexData(i);
      V.connectEdge(this->id);
    }
  };
  /**
   * @brief Get a pointer to the local vertex local_number of the quadrilateral
   * face.
   *
   * @param pointers object containing the pointers to global data.
   * @param local_number local number of the vertex.
   * @return Geometry::Vertex* pointer to the vertex.
   */
  auto getVertex(indexType local_number)
      -> Geometry::VertexData * override {
    return this->m_verts_pointers[local_number];
  };
  /**
   * @brief Get the global vertex numbers of the geometric element.
   *
   * @param vertsOut[out], a vector with the global numbers of the vertices of
   * type indexType and of size 8.
   */
  auto getVertexNumbers() -> std::vector<indexType> override {
    std::vector<indexType> nums(this->m_verts.begin(), this->m_verts.end());
    return nums;
  };
  /** @brief Get the global number of the vertex "number" of the object.
   *
   * @return indexType Number of Vertex.
   */
  auto getVertexNumber(indexType number) -> indexType override {
    return this->m_verts[number];
  };
  /** @brief Get an array of VertexData pointer to the vertices.
   *
   * @return indexType Number of Vertex.
   */
  void getVerts(std::vector<VertexData *> &vertsOut) override {
    vertsOut = std::vector<VertexData *>(this->m_verts_pointers.begin(),
                                               this->m_verts_pointers.end());
  };
  /**
   * @brief Get the Number Of Vertices of the geometric object, return
   * m_numberOfVerts
   *
   * @return indexType, returns m_numberOfVerts.
   */
  auto getNumberOfVerts() -> indexType override { return this->m_numberOfVerts; };


  /**
   * @brief Updates the pointers to the geometry objects
   *
   */
  void set_geometry_pointers(GeometryData &geoData) override {
    for (auto i = 0; i < this->m_numberOfVerts; ++i) {
      if (this->m_verts[i] != -1)
        this->m_verts_pointers[i] = &geoData.getVertexData(this->m_verts[i]);
    }
  };

  void print(spdlog::logger &Logger) override {

    Logger.debug("{} id: {:8d}", self::getName(), this->id);
    Logger.debug("Vertices: {}", fmt::join(this->m_verts, " "));

    static_cast<self *>(this)->printEqInfo(Logger);
  };

  /** @brief Check if the face has the three vertices.
   *
   * @param[in] v1 Global number of a corner vertex.
   * @param[in] v2 Global number of a corner vertex.
   * @return true If the Edge has the two vertices.
   * @return false If the Edge does not have the two vertices.
   */
  auto hasVertices(indexType v1, indexType v2) -> bool override {
    return (this->m_verts[0] == v1 && this->m_verts[this->m_numberOfVerts-1] == v2) ||
           (this->m_verts[0] == v2 && this->m_verts[this->m_numberOfVerts-1] == v1);
  };

protected:
};

} // namespace HierAMuS::Geometry
