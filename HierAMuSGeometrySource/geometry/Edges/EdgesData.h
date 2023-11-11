// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MeshIdNodeList.h"

#include "datatypes.h"
#include "geometry/GeometryBaseData.h"
#include "geometry/GeometryShape.h"

#include <geometry/GeometryTypes.h>
#include <types/MatrixTypes.h>

#include <vector>

namespace HierAMuS {
class vtkPlotInterface;
}

namespace HierAMuS::Geometry {
class EdgesRuntime;

class VertexData;
class EdgesData : public GeometryBaseData {
public:
  EdgesData();
  ~EdgesData() override;

  virtual auto getRuntimeObject(GeometryData &geoData)
      -> std::shared_ptr<EdgesRuntime> = 0;

  auto getGroupType() -> const GeometryTypes & override;
  virtual auto getEdgeOrientation(indexType startNode, indexType endNode)
      -> prec = 0;

  virtual auto getA1Vector(IntegrationPoint &integration_point)
      -> Types::Vector3<prec> {
    throw std::runtime_error(
        "Function getA1Vector not implemented for the used Edge object!");
    return {};
  };

  virtual auto getJacobian(prec xi) -> prec {
    return prec(0);
  };

  virtual auto getJacobian(IntegrationPoint &IntegrationPt)
      -> prec {
    throw std::runtime_error(
        "Function getJacobian not implemented for the used Edge object!");
    return {};
  }

  /**
   * @brief Get the physical Coordinates at point xi of parametric coordinate
   * system.
   *
   * @param[in] pointers Global data collection.
   * @param[in] xi Parameteric position between -1 and 1.
   * @return Types::Vector3<prec> Physical coordinates x y z.
   */
  virtual auto getCoordinates(prec xi)
      -> Types::Vector3<prec> = 0;
  virtual auto getCoordinates(IntegrationPoint &IntPoint)
      -> Types::Vector3<prec> = 0;

  /**
   * @brief Get the Vertex "number" of object.
   *
   * @param[in] pointers Global data collection.
   * @param[in] number Local Vertex number.
   * @return Geometry::Vertex* Pointer to vertex object.
   */
  virtual auto getVertex(indexType number)
      -> VertexData * = 0;

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
  /**
   * @brief Get the global vertex numbers of the geometric element.
   *
   * @param vertsOut[out], a vector with the global numbers of the vertices of
   * type indexType and of size 8.
   */
  virtual auto getVertexNumbers() -> std::vector<indexType> = 0;

  // H1 Shapes
  virtual void setH1Shapes(indexType meshId,
                           indexType order, NodeTypes type) = 0;
  virtual void setH1ShapesInternal(indexType meshId, indexType order,
                                   NodeTypes type) = 0;

  virtual void getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) = 0;

  virtual void getH1DofsInternal(std::vector<DegreeOfFreedom *> &Dofs,
                                 indexType meshID, indexType order) = 0;

  virtual auto getH1NodesList(indexType meshID,
                              indexType order) -> MeshIdNodeList = 0;

  virtual auto getH1Nodes(indexType meshID,
                          indexType order) -> std::vector<GenericNodes *> = 0;

  virtual auto getH1NodesInternal(indexType meshID, indexType order)
      -> std::vector<GenericNodes *> = 0;

  virtual void getH1Shapes(indexType order,
                           Types::VectorX<prec> &shape,
                           Types::VectorX<prec> &shapeDerivative, prec xsi) {
    EdgesData::throwError(
        "Error when calling getH1Shapes for Edges, not implemented!");
  };

  virtual void getH1ShapesInternal(indexType order,
                                   Types::VectorX<prec> &shape,
                                   Types::VectorX<prec> &shapeDerivative,
                                   prec xsi) {
    EdgesData::throwError(
        "Error when calling getH1ShapesInternal for Edges, not implemented!");
  };

  virtual auto getH1Shapes(indexType order,
                           IntegrationPoint &integration_point) -> H1Shapes = 0;

  virtual auto getH1ShapesInternal(indexType order,
                                   IntegrationPoint &integration_point)
      -> H1Shapes {
    EdgesData::throwError(
        "Error when calling getH1Shapes for Edges, not implemented!");
    return H1Shapes(0, 0);
  };

  // HDivShapes
  void setHDivShapes(indexType meshid,
                     indexType order, NodeTypes type) {
    EdgesData::throwError(
        "Error when calling setHDivShapes for Faces, not implemented!");
  };

  void getHDivDofs(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                   indexType order, NodeTypes type) {
    EdgesData::throwError(
        "Error when calling getHDivDofs for Faces, not implemented!");
  };
  void getHDivShapes(indexType order,
                     Types::Matrix22<prec> &jacobi, Types::VectorX<prec> &shape,
                     Types::Matrix2X<prec> &dshape, prec xi, prec eta) {
    HierAMuS::Geometry::EdgesData::throwError(
        "Error when calling getHDivShapes for Faces, not implemented!");
  };

  virtual auto getDirectionVector()
      -> Types::Vector3<prec> {
    return Types::Vector3<prec>(0);
  };

  void setAllNodeBoundaryConditionMeshId(indexType meshId,
                                         indexType dof) override;

  /**
   * @brief Changes orientation of the face (reverts the face normal vector).
   */
  virtual void flip() = 0;

  

  virtual void set_geometry_pointers(GeometryData &geoData) = 0;

  virtual void getVerts(std::vector<VertexData *> &vertsOut) = 0;

  virtual void setVerts(
      GeometryData &geoData,
      std::vector<indexType> &vertsIn) = 0;
  virtual auto getNumberOfVerts() -> indexType = 0;


private:
  static const HierAMuS::Geometry::GeometryTypes type;
  static void throwError(const std::string &msg) {
    throw std::runtime_error(msg);
  };
};

} // namespace HierAMuS::Geometry
