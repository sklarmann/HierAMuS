// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "datatypes.h"
#include "MatrixTypes.h"
#include "geometry/GeometryTypes.h"
#include "geometry/VertexRuntime.h"
#include <stack>
#include <vector>
#include <unordered_map>
#include <memory>
#include "spdlog/spdlog.h"


//#include "geometry/VertexRuntime.h"

namespace HierAMuS {
class EquationHandler;
}

namespace HierAMuS::Geometry {
class GeometryBaseData;
class VertexData;
class VertexRuntime;
class EdgesData;
class EdgesRuntime;
class FacesData;
class FacesRuntime;
class VolumesData;
class VolumesRuntime;
class Special;
template <class VertexData> class GeometryListSingle;
template <class EdgesData> class GeometryList;
template <class FacesData> class GeometryList;
template <class VolumesData> class GeometryList;
class GeometryData {
public:
  GeometryData();
  virtual ~GeometryData();
  auto requestNewVert() -> indexType;
  

  auto getVertexData(indexType num) -> VertexData &;
  auto getVertexRuntime(indexType num) -> VertexRuntime&;
  auto getVertexRuntimePointer(indexType num)
      -> VertexRuntime*;
  auto getEdgeData(indexType num) -> EdgesData &;
  auto getEdgeRuntime(indexType num) -> std::shared_ptr<EdgesRuntime>;
  auto getFaceData(indexType num) -> FacesData *;
  auto getFaceRuntime(indexType num) -> std::shared_ptr<FacesRuntime>;
  auto getVolumeData(indexType num) -> VolumesData *;
  auto getVolumeRuntime(indexType num)
      -> std::shared_ptr<VolumesRuntime>;
  auto getSpecial(indexType num) -> Special *;

  auto getVertexClosestTo(Types::Vector3<prec> point) -> VertexData &;

  auto requestNewGeometryObject(EquationHandler &eqHandler, GeometryTypes type)
      -> indexType;
  auto requestNewGeometryObject(EquationHandler &eqHandler, GeometryTypes type, indexType number)
      -> indexType;
  auto getGeometryElement(GeometryTypes type, indexType num) -> GeometryBaseData *;

  auto getNumberOfVertices() -> indexType;
  auto getNumberOfEdges() -> indexType;

  auto getLastVertexNumber() -> indexType;
  auto getNextVertexNumber() -> indexType;

  void getGeometricElementInPlane(std::vector<prec> normal,
                                  std::vector<prec> point, GeometryTypes type,
                                  std::stack<GeometryBaseData *> &elems);

  auto getVerticesInPlane(const Types::Vector3<prec> &normal,
                          const Types::Vector3<prec> &point)
      -> std::vector<Geometry::VertexData *>;
  auto getEdgesInPlane(const Types::Vector3<prec> &normal,
                       const Types::Vector3<prec> &point)
      -> std::vector<Geometry::EdgesData *>;
  auto getFacesInPlane(const Types::Vector3<prec> &normal,
                       const Types::Vector3<prec> &point)
      -> std::vector<Geometry::FacesData *>;

  void print(spdlog::logger &Logger);

  void checkUpdate(EquationHandler &eqHandler);

  auto getEdgeNumberByVerts(indexType vert1, indexType vert2) -> indexType;
  auto getFaceNumberByVerts(indexType vert1, indexType vert2, indexType vert3)
      -> indexType;

  auto getxMax() -> Types::Vector3<prec>;
  auto getxMin() -> Types::Vector3<prec>;

  /** @brief Sorts the slave face numbers such that they have the same local x1,
   * x2 coordinates and are just shifted by x3. Reorients the slave faces such
   * that their local coordinate system matches the one of the master faces,
   * this is required for periodic boundary conditions with hierarchical higher
   * order shape functions.
   *
   *
   * @param[inout] masterFaces vector with the master face numbers.
   * @param[inout] slaveFaces vector with the slave face numbers.
   */
  void sortReorientFacesPeriodicBC(std::vector<indexType> &masterFaces,
                                   std::vector<indexType> &slaveFaces);

  void createRuntimeObjects();
  void updateRuntimeObjectEquations();

private:
  std::shared_ptr<GeometryListSingle<VertexData>> m_vertices;
  std::shared_ptr<GeometryList<EdgesData>> m_edges;
  std::shared_ptr<GeometryList<FacesData>> m_faces;
  std::shared_ptr<GeometryList<VolumesData>> m_volumes;
  std::shared_ptr<GeometryList<Special>> m_special;

  //std::shared_ptr<GeometryListSingle<VertexRuntime>> m_vertices_runtime;
  std::map<indexType, std::shared_ptr<VertexRuntime>> m_vertices_runtime;
  std::map<indexType, std::shared_ptr<EdgesRuntime>> m_edges_runtime;
  std::map<indexType, std::shared_ptr<FacesRuntime>> m_faces_runtime;
  std::map<indexType, std::shared_ptr<VolumesRuntime>> m_volumes_runtime;

  Types::Vector3<prec> xMin, xMax;

};

} // namespace HierAMuS::Geometry
