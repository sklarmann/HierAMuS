// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "geometry/Faces/FacesData.h"

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

template <bool, indexType numVertices, indexType numEdges>
class FacesDataInterfaceVariables : public FacesData {};

// Static sized storage
template <indexType numVertices, indexType numEdges>
class FacesDataInterfaceVariables<true, numVertices, numEdges>
    : public FacesData {
public:
  FacesDataInterfaceVariables()
      : FacesData(), m_verts(initArray<numVertices, indexType, -1>()),
        m_edges(initArray<numEdges, indexType, -1>()){};

protected:
  inline static const indexType m_numberOfVerts = numVertices;
  inline static const indexType m_numberOfEdges = numEdges;

  std::array<indexType, numVertices> m_verts;
  std::array<indexType, numEdges> m_edges;

  std::array<VertexData *, numVertices> m_verts_pointers;
  std::array<EdgesData *, numEdges> m_edges_pointers;
};

// Dynamic sized storage
template <indexType numVertices, indexType numEdges>
class FacesDataInterfaceVariables<false, numVertices, numEdges>
    : public FacesData {
public:
  FacesDataInterfaceVariables()
      : FacesData(), m_numberOfVerts(0), m_numberOfEdges(0){};

  void setNumberOfVertices(indexType numberOfVertices) {
    m_numberOfVerts = numberOfVertices;
    m_verts.resize(numberOfVertices);
    m_verts_pointers.resize(numberOfVertices);
    std::fill(m_verts.begin(), m_verts.end(), -1);
  }
  void setNumberOfEdges(indexType numberOfEdges) {
    m_numberOfEdges = numberOfEdges;
    m_edges.resize(numberOfEdges);
    m_edges_pointers.resize(numberOfEdges);
    std::fill(m_edges.begin(), m_edges.end(), -1);
  }

protected:
  indexType m_numberOfVerts;
  indexType m_numberOfEdges;

  std::vector<indexType> m_verts;
  std::vector<indexType> m_edges;

  std::vector<VertexData *> m_verts_pointers;
  std::vector<EdgesData *> m_edges_pointers;
};

template <indexType numVertices, indexType numEdges, class self>
class FacesDataInterface
    : public std::conditional_t<
          numVertices <= 0 || numEdges <= 0,
          FacesDataInterfaceVariables<false, numVertices, numEdges>,
          FacesDataInterfaceVariables<true, numVertices, numEdges>>
//  std::conditional_t<
//          numVertices <= 0 || numEdges <= 0,
//          FacesDataInterfaceVariables<false, numVertices, numEdges>,
//          FacesDataInterfaceVariables<true, numVertices, numEdges>>
{
  using baseBool = typename std::conditional_t<
          numVertices <= 0 || numEdges <= 0,
          std::integral_constant<bool, false>,
          std::integral_constant<bool, true>>;
public:
  FacesDataInterface() : FacesDataInterfaceVariables<baseBool::value, numVertices, numEdges>(){};
  ~FacesDataInterface() = default;
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
      V.connectFace(this->id);
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
   * @brief Set the Edge numbers of the geometric element.
   * They need to be given according to the numbering in the picture above.
   *
   * @param vertsIn[in], a vector with the global numbers of the Edges of type
   * indexType and of size 12.
   */
  void setEdges(const std::vector<indexType> &edgesIn) override {
    std::copy_n(edgesIn.begin(), this->m_numberOfEdges, this->m_edges.begin());
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
  /** @brief Get an array of GeometryBaseData pointer to the vertices.
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

  /** @brief Get the global number of the edge "number" of the object.
   *
   * @return indexType Number of edge.
   */
  auto getEdgeNumber(indexType localNumber) -> indexType override{
    return this->m_edges[localNumber];
  };
  /** @brief Get the global number of the edge "number" of the object.
   *
   * @return indexType Number of edge.
   */
  auto getEdge(indexType localNumber) -> EdgesData * override {
    return this->m_edges_pointers[localNumber];
  };
  /**
   * @brief Get the global edge numbers of the geometric element.
   *
   * @param edgesOut[out], a vector with the global numbers of the Edges of type
   * indexType and of size 12.
   */
  void getEdgeNumbers(std::vector<indexType> &edgesOut) override {
    edgesOut = std::vector<indexType>(this->m_edges.begin(), this->m_edges.end());
  };
  /**
   * @brief Get the global edge numbers of the object.
   *
   * @param[out] edgesOut, vector of 3 integers of type indexType with global
   * edge numbers of the element.
   */
  auto getEdgeNumbers() -> std::vector<indexType> override {
    return std::vector<indexType>(this->m_edges.begin(), this->m_edges.end());
  };
  /**
   * @brief Get a vector with pointers of type Geometry::Base to the edges of
   * the geometric element.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param edgesOut[out], a vector of size 12 with pointers to the edges of the
   * geometric element.
   */
  void getEdges(std::vector<EdgesData *> &edgesOut) override {
    edgesOut = std::vector<EdgesData *>(this->m_edges_pointers.begin(),
                                               this->m_edges_pointers.end());
  };
  /**
   * @brief Get the Number of Edges of the element, returns 12.
   *
   * @return indexType, returns 12.
   */
  auto getNumberOfEdges() -> indexType override { return this->m_numberOfEdges; };

  /**
   * @brief Updates the pointers to the geometry objects
   *
   */
  void set_geometry_pointers(GeometryData &geoData) override {
    for (auto i = 0; i < this->m_numberOfVerts; ++i) {
      if (this->m_verts[i] != -1)
        this->m_verts_pointers[i] = &geoData.getVertexData(this->m_verts[i]);
    }
    for (auto i = 0; i < this->m_numberOfEdges; ++i) {
      if (this->m_edges[i] != -1)
        this->m_edges_pointers[i] = &geoData.getEdgeData(this->m_edges[i]);
    }
  };

  void print(spdlog::logger &Logger) override {

    Logger.debug("{} id: {:8d}", self::getName(), this->id);
    Logger.debug("Vertices: {}", fmt::join(this->m_verts, " "));
    Logger.debug("Edges:    {}", fmt::join(this->m_edges, " "));

    static_cast<self *>(this)->printEqInfo(Logger);
  };

  /** @brief Check if the face has the three vertices.
   *
   * @param[in] v1 Global number of a corner vertex.
   * @param[in] v2 Global number of a corner vertex.
   * @param[in] v3 Global number of a corner vertex.
   * @return true If the face has the three vertices.
   * @return false If the face does not have the three vertices.
   */
  auto hasVertices(indexType v1, indexType v2, indexType v3) -> bool override {
    if (std::find(this->m_verts.begin(), this->m_verts.end(), v1) != this->m_verts.end()) {
      if (std::find(this->m_verts.begin(), this->m_verts.end(), v2) != this->m_verts.end()) {
        if (std::find(this->m_verts.begin(), this->m_verts.end(), v3) != this->m_verts.end()) {
          return true;
        }
      }
    }
    return false;
  };

  void getNodes(std::vector<GenericNodes *> &nodeVector,
                                indexType meshId) override {
    nodeVector.clear();
    for (auto i : this->m_verts_pointers) {
      auto tempNodes = i->getNodesOfSet(meshId);
      nodeVector.insert(nodeVector.end(), tempNodes.begin(), tempNodes.end());
    }
    for (auto i : this->m_edges_pointers) {
      auto tempNodes = i->getNodesOfSet(meshId);
      nodeVector.insert(nodeVector.end(), tempNodes.begin(), tempNodes.end());
    }

    auto tempNodes = this->getNodesOfSet(meshId);
    nodeVector.insert(nodeVector.end(), tempNodes.begin(), tempNodes.end());


  };

protected:
};

} // namespace HierAMuS::Geometry
