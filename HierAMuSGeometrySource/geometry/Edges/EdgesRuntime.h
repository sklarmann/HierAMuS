// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MeshIdNodeList.h"

#include "datatypes.h"
#include "geometry/GeometryBaseData.h"
#include "geometry/GeometryShape.h"
#include <geometry/GeometryBaseRuntime.h>

#include <geometry/GeometryTypes.h>
#include <types/MatrixTypes.h>

#include <vector>

namespace HierAMuS {
class vtkPlotInterface;
} // namespace HierAMuS

namespace HierAMuS::Geometry {
class VertexRuntime;
class EdgesData;
class EdgeH1ShapesInterface;

class VertexData;
class EdgesRuntime : public GeometryBaseRuntime {
public:
  EdgesRuntime(EdgesData &data_element);
  ~EdgesRuntime() override;
  auto getGroupType() -> const GeometryTypes & override;
  virtual auto getEdgeOrientation(indexType startNode, indexType endNode)
      -> prec = 0;

  virtual auto getH1Edge() -> EdgeH1ShapesInterface * { return NULL; };

  virtual auto getA1Vector(IntegrationPoint &integration_point)
      -> Types::Vector3<prec> {
    throw std::runtime_error(
        "Function getA1Vector not implemented for the used Edge object!");
    return {};
  };

  virtual auto getJacobian(IntegrationPoint &ip) -> prec = 0;

  /**
   * @brief Get the physical Coordinates at point xi of parametric coordinate
   * system.
   *
   * @param[in] pointers Global data collection.
   * @param[in] xi Parametric position between -1 and 1.
   * @return Types::Vector3<prec> Physical coordinates x y z.
   */
  virtual auto getCoordinates(prec xi) -> Types::Vector3<prec> = 0;
  virtual auto getCoordinates(IntegrationPoint &IntPoint)
      -> Types::Vector3<prec> = 0;

  /**
   * @brief Get the Vertex "number" of object.
   *
   * @param[in] pointers Global data collection.
   * @param[in] number Local Vertex number.
   * @return Geometry::Vertex* Pointer to vertex object.
   */
  virtual auto getVertex(indexType number) -> VertexRuntime & = 0;

  virtual void getVertices(std::vector<VertexRuntime *> &verticesOut) = 0;

  /** @brief Check if the edge has the two vertices.
   *
   * @param[in] v1 Global number of the start vertex.
   * @param[in] v2 Global number of the end vertex.
   * @return true If the edge has the two vertices.
   * @return false If the edge does not have the two vertices.
   */
  virtual auto hasVertices(indexType v1, indexType v2) -> bool = 0;

  /** @brief Get the global number of the vertex "number" of the object.
   *
   * @return indexType Number of Vertex.
   */
  virtual auto getVertexNumber(indexType number) -> indexType = 0;

  // HDivShapes
  void setHDivShapes(indexType meshid, indexType order, NodeTypes type) {
    EdgesRuntime::throwError(
        "Error when calling setHDivShapes for Faces, not implemented!");
  };

  void getHDivDofs(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                   indexType order, NodeTypes type) {
    EdgesRuntime::throwError(
        "Error when calling getHDivDofs for Faces, not implemented!");
  };
  void getHDivShapes(indexType order, Types::Matrix22<prec> &jacobi,
                     Types::VectorX<prec> &shape, Types::Matrix2X<prec> &dshape,
                     prec xi, prec eta) {
    HierAMuS::Geometry::EdgesRuntime::throwError(
        "Error when calling getHDivShapes for Faces, not implemented!");
  };

  virtual auto getDirectionVector() -> Types::Vector3<prec> {
    return Types::Vector3<prec>(0);
  };

  void setAllNodeBoundaryConditionMeshId(indexType meshId,
                                         indexType dof) override;

  /**
   * @brief Changes orientation of the face (reverts the face normal vector).
   */
  virtual void flip() = 0;

  // Paraview Part
  virtual void geometryToParaview(vtkPlotInterface &paraviewAdapter,
                                  indexType mainMesh, indexType subMesh){};
  virtual void computeWeightsParaview(vtkPlotInterface &paraviewAdapter,
                                      indexType mainMesh, indexType subMesh){};
  virtual void H1SolutionToParaview(vtkPlotInterface &paraviewAdapter,
                                    indexType mainMesh, indexType subMesh,
                                    indexType order,
                                    Types::VectorX<prec> solution,
                                    std::string &name){};
  virtual void H1DataToParaview(vtkPlotInterface &paraviewAdapter,
                                indexType mainMesh, indexType subMesh,
                                Types::VectorX<prec> &Data,
                                indexType numberComponents, indexType order,
                                std::string &name){};
  virtual void projectDataToParaviewVertices(
      vtkPlotInterface &paraviewAdapter, indexType mainMesh, indexType subMesh,
      indexType order, IntegrationPoint &IntegrationPt,
      Types::VectorX<prec> &data, indexType numberComponents,
      std::string name){}; // move into specialized classes

  virtual void set_geometry_pointers(GeometryData &geoData) = 0;

  // Getters and setters
  virtual auto getNumberOfVerts() -> indexType = 0;

  virtual auto getVertexNumbers() -> std::vector<indexType> = 0;

private:
  static const HierAMuS::Geometry::GeometryTypes type;
  static void throwError(const std::string &msg) {
    throw std::runtime_error(msg);
  };
};

} // namespace HierAMuS::Geometry
