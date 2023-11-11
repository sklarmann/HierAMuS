// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <geometry/Volumes/VolumesData.h>

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
class VolumesData;
class FacesRuntime;

template <indexType numVertices, indexType numEdges, indexType numFaces,
          class self>
class VolumesDataInterface : public VolumesData {
public:
  VolumesDataInterface()
      : VolumesData(), m_verts(initArray<m_numberOfVerts, indexType, -1>()),
        m_edges(initArray<m_numberOfEdges, indexType, -1>()),
        m_faces(initArray<m_numberOfFaces, indexType, -1>()){};
  ~VolumesDataInterface() = default;
  /**
   * @brief Set the Vertex numbers of the geometric element.
   * They need to be given according to the numbering in the picture above.
   *
   * @param vertsIn[in], a vector with the global numbers of the vertices of
   * type indexType and of size 8.
   */
  void setVerts(GeometryData &geoData,
                std::vector<indexType> &vertsIn) override {
    std::copy_n(vertsIn.begin(), m_numberOfVerts, m_verts.begin());
  };
  /**
   * @brief Set the Edge numbers of the geometric element.
   * They need to be given according to the numbering in the picture above.
   *
   * @param vertsIn[in], a vector with the global numbers of the Edges of type
   * indexType and of size 12.
   */
  void setEdges(const std::vector<indexType> &edgesIn) override {
    std::copy_n(edgesIn.begin(), m_numberOfEdges, m_edges.begin());
  };
  /**
   * @brief Set the Face numbers of the geometric element.
   * They need to be given according to the numbering in the picture above.
   *
   * @param vertsIn[in], a vector with the global numbers of the Faces of type
   * indexType and of size 6.
   */
  void setFaces(std::vector<indexType> &facesIn) override {
    std::copy_n(facesIn.begin(), m_numberOfFaces, m_faces.begin());
  };
  /**
   * @brief Get the global vertex numbers of the geometric element.
   *
   * @param vertsOut[out], a vector with the global numbers of the vertices of
   * type indexType and of size 8.
   */
  auto getVertexNumbers() -> std::vector<indexType> override {
    std::vector<indexType> nums(m_verts.begin(), m_verts.end());
    return nums;
  };
  /** @brief Get the global number of the vertex "number" of the object.
   *
   * @return indexType Number of Vertex.
   */
  auto getVertexNumber(indexType number) -> indexType override {
    return m_verts[number];
  };
  /** @brief Get an array of GeometryBaseData pointer to the vertices.
   *
   * @return indexType Number of Vertex.
   */
  void getVerts(std::vector<VertexData *> &vertsOut) override {
    vertsOut = std::vector<VertexData *>(m_verts_pointers.begin(),
                                         m_verts_pointers.end());
  };
  /**
   * @brief Get the Number Of Vertices of the geometric object, return
   * m_numberOfVerts
   *
   * @return indexType, returns m_numberOfVerts.
   */
  auto getNumberOfVerts() -> indexType override { return m_numberOfVerts; };

  /** @brief Get the global number of the edge "number" of the object.
   *
   * @return indexType Number of edge.
   */
  auto getEdgeNumber(indexType localNumber) -> indexType override {
    return m_edges[localNumber];
  };
  /** @brief Get the global number of the edge "number" of the object.
   *
   * @return indexType Number of edge.
   */
  auto getEdge(indexType localNumber) -> EdgesData * override {
    return m_edges_pointers[localNumber];
  };
  /**
   * @brief Get the global edge numbers of the geometric element.
   *
   * @param edgesOut[out], a vector with the global numbers of the Edges of type
   * indexType and of size 12.
   */
  void getEdgeNumbers(std::vector<indexType> &edgesOut) override {
    edgesOut = std::vector<indexType>(m_edges.begin(), m_edges.end());
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
    edgesOut = std::vector<EdgesData *>(m_edges_pointers.begin(),
                                               m_edges_pointers.end());
  };
  /**
   * @brief Get the Number of Edges of the element, returns 12.
   *
   * @return indexType, returns 12.
   */
  auto getNumberOfEdges() -> indexType override { return m_numberOfEdges; };

  /** @brief Get the global number of the face "number" of the object.
   *
   * @return indexType Number of face.
   */
  auto getFaceNumber(indexType number) -> indexType override { return m_faces[number]; };
  /**
   * @brief Get the global face numbers of the geometric element.
   *
   * @param facesOut[out], a vector with the global numbers of the Faces of type
   * indexType and of size 6.
   */
  void getFaceNumbers(std::vector<indexType> &facesOut) override {
    facesOut = std::vector<indexType>(m_faces.begin(), m_faces.end());
  };

  /**
   * @brief Get a vector with pointers of type Geometry::Base to the faces of
   * the geometric element.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param facesOut[out], a vector of size 6 with pointers to the faces of the
   * geometric element.
   */
  void getFaces(std::vector<GeometryBaseData *> &faces) override {
    faces = std::vector<GeometryBaseData *>(m_faces_pointers.begin(),
                                            m_faces_pointers.end());
  };
  /**
   * @brief Get the Number Of Faces object
   *
   * @return indexType Number of faces as an indexType
   */
  auto getNumberOfFaces() -> indexType override { return m_numberOfFaces; };

  /**
   * @brief Updates the pointers to the geometry objects
   *
   */
  void set_geometry_pointers(GeometryData &geoData) override {
    for (auto i = 0; i < m_numberOfVerts; ++i) {
      if (m_verts[i] != -1)
        m_verts_pointers[i] = &geoData.getVertexData(m_verts[i]);
    }
    for (auto i = 0; i < m_numberOfEdges; ++i) {
      if (m_edges[i] != -1)
        m_edges_pointers[i] = &geoData.getEdgeData(m_edges[i]);
    }
    for (auto i = 0; i < m_numberOfFaces; ++i) {
      if (m_faces[i] != -1)
        m_faces_pointers[i] = geoData.getFaceData(m_faces[i]);
    }
  };

  void print(spdlog::logger &Logger) override {

    Logger.debug("{} id: {:8d}", self::getName(), this->id);
    Logger.debug("Vertices: {}", fmt::join(m_verts, " "));
    Logger.debug("Edges:    {}", fmt::join(m_edges, " "));
    Logger.debug("Faces:    {}", fmt::join(m_faces, " "));

    static_cast<self *>(this)->printEqInfo(Logger);
  };


  void getNodes(std::vector<GenericNodes *> &nodeVector, indexType meshId) override {
    nodeVector.clear();
    for (auto i : m_verts_pointers) {
      auto tempNodes = i->getNodesOfSet(meshId);
      nodeVector.insert(nodeVector.end(), tempNodes.begin(), tempNodes.end());
    }
    for (auto i : m_edges_pointers) {
      auto tempNodes = i->getNodesOfSet(meshId);
      nodeVector.insert(nodeVector.end(), tempNodes.begin(), tempNodes.end());
    }
    for (auto i : m_faces_pointers) {
      auto tempNodes = i->getNodesOfSet(meshId);
      nodeVector.insert(nodeVector.end(), tempNodes.begin(), tempNodes.end());
    }

    auto tempNodes = this->getNodesOfSet(meshId);
    nodeVector.insert(nodeVector.end(), tempNodes.begin(), tempNodes.end());
  };

protected:
  static constexpr indexType m_numberOfVerts = numVertices;
  static constexpr indexType m_numberOfEdges = numEdges;
  static constexpr indexType m_numberOfFaces = numFaces;

  std::array<indexType, m_numberOfVerts> m_verts;
  std::array<indexType, m_numberOfEdges> m_edges;
  std::array<indexType, m_numberOfFaces> m_faces;

  // Erster Schritt fuer spaetere Separierung in Index- und Berechnungsklasse
  std::array<VertexData *, m_numberOfVerts> m_verts_pointers;
  std::array<EdgesData *, m_numberOfEdges> m_edges_pointers;
  std::array<FacesData *, m_numberOfFaces> m_faces_pointers;
};

} // namespace HierAMuS::Geometry
