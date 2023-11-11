// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <geometry/GeometryBaseData.h>

#include "geometry/GeometryShape.h"
#include <geometry/GeometryTypes.h>
#include <types/MatrixTypes.h>

#include <vector>

namespace HierAMuS {
class EquationHandler;
}

namespace HierAMuS::Geometry {
class VolumesRuntime;
class GeometryData;
class VertexData;
class EdgesData;
class VolumesData : public GeometryBaseData {
public:
  VolumesData();
  ~VolumesData() override;
  auto getGroupType() -> const GeometryTypes & override;

  virtual auto getRuntimeObject(GeometryData &geoData)
      -> std::shared_ptr<VolumesRuntime> = 0;
  /**
   * @brief Returns the number of edges of the element.
   * @return Number of edges as an indexType
   */
  virtual auto getNumberOfEdges() -> indexType = 0;
  /**
   * @brief Get the Number Of Faces object
   *
   * @return indexType Number of faces as an indexType
   */
  virtual auto getNumberOfFaces() -> indexType = 0;

  virtual void setFaces(std::vector<indexType> &facesIn) = 0;

  // H1 Shapes
  virtual void setH1Shapes(indexType meshId,
                           indexType order, NodeTypes type) = 0;
  virtual void setH1ShapesInternal(indexType meshId, indexType order,
                                   NodeTypes type) = 0;

  // Checking functions
  /**
   * @brief Checks if the element is completely defined.
   *
   * Checks if the element is completely defined. If not, it will search or
   * create the necessary additional geometry elements and add them to the
   * element.
   *
   * @param geoData Pointer to the geometry data object
   */
  virtual void checkUpdateElement(EquationHandler &eqHandler, GeometryData &geoData) = 0;

  virtual void set_geometry_pointers(GeometryData &geoData) = 0;

  virtual auto getVertexNumbers() -> std::vector<indexType> = 0;
  /** @brief Get the global number of the vertex "number" of the object.
   *
   * @return indexType Number of Vertex.
   */
  virtual auto getVertexNumber(indexType number) -> indexType = 0;
  virtual void getVerts(std::vector<VertexData *> &vertsOut) = 0;

  virtual auto getEdgeNumber(indexType localNumber) -> indexType = 0;

  virtual void getEdgeNumbers(std::vector<indexType> &edgesOut) = 0;
  virtual void getEdges(std::vector<EdgesData *> &edgesOut) = 0;
  virtual auto getEdge(indexType localNumber) -> EdgesData * = 0;

  virtual auto getFaceNumber(indexType localNumber) -> indexType = 0;

  virtual void getFaceNumbers(std::vector<indexType> &facesOut) = 0;

  virtual void getFaces(std::vector<GeometryBaseData *> &faces) = 0;

  virtual void setEdges(const std::vector<indexType> &edgesIn) = 0;

  virtual void setVerts(GeometryData &geoData, std::vector<indexType> &vertsIn) = 0;
  virtual auto getNumberOfVerts() -> indexType = 0;

private:
  static const GeometryTypes type;
  static void throwError(const std::string &msg) {
    throw std::runtime_error(msg);
  };
};

} // namespace HierAMuS::Geometry
