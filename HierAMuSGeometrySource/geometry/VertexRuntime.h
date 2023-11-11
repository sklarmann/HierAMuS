// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "datatypes.h"
#include <array>

#include "geometry/GeometryBaseData.h"
#include "geometry/GeometryBaseRuntime.h"
#include <geometry/GeometryTypes.h>


namespace HierAMuS::Geometry {
class VertexData;

class VertexRuntime : public GeometryBaseRuntime {
public:
  VertexRuntime(VertexData &data);
  VertexRuntime(VertexRuntime &other);
  VertexRuntime(VertexRuntime &&other);
  VertexRuntime &operator=(VertexRuntime &other);
  VertexRuntime &operator=(VertexRuntime &&other);
  ~VertexRuntime() override;
  /**
   * @brief Routine to set the coordinates of the vertex.
   *
   * @param x x Coordinate
   * @param y y Coordinate
   * @param z z Coordinate
   */
  void setCoordinates(prec x = 0, prec y = 0, prec z = 0);
  /**
   * @brief Routine to set coordinates of the vertex.
   *
   * @param[in] coorIn Vector3 with the x, y, z coordinates of the vertex.
   */
  void setCoordinates(const Types::Vector3<prec> &coorIn);
  void print(spdlog::logger &Log) override;
  auto getType() -> const HierAMuS::Geometry::GeometryTypes & override;
  auto getCoordinates() -> Types::Vector3<prec>;

  /**
   * @brief Set Boundary conditions on the element
   *
   * @param[in] meshId Mesh id of the nodes.
   * @param[in] order Order of the shape function.
   * @param[in] shapeType Type of the shape function to restrict.
   * @param[in] dofs Degree of freedom list to set boundary conditions.
   * @param[in] set If true, boundary conditions are overridden, otherwise only 1
   * is considered.
   */
  void setBoundaryCondition(indexType meshId,
                            indexType order, ShapeFunctionTypes shapeType,
                            Types::Vector3<indexType> &dofs, bool set) override;

  /**
   * @brief Sets load using a Geometric Object.
   *
   * @param[in] meshid Mesh id of the degrees of freedom to add the load.
   * @param[in] shapeType Shape function Type (H0 or H1).
   * @param[in] shapeOrder Order of the shape functions to approximate delta u.
   * @param[in] Loads Vector with loads, length must be a multiple of 3,
   * depending on the length, the load interpolation function is chosen.
   * @param[in] propNumber Proportional load function number for the load.
   * @param[in] direction If local is a function, this vector estimates if the
   * integration is in positive or negative direction of an edge.
   * @param[in] local If true, the loads are treated in the local coordinate
   * system.
   * @param[in] add If true, the loads will be added to the current loads
   * otherwise loads will be overridden.
   */
  void setLoad(LoadList &loadlist, indexType meshid,
               ShapeFunctionTypes shapeType, indexType shapeOrder,
               Types::VectorX<prec> &Loads, indexType propNumber,
               Types::VectorX<prec> &direction, bool local, bool add) override;

  /**
   * @brief Sets the solution using a Geometric Object.
   * If boundary conditions are not set, they will be set, if the prescribed
   * value is not zero.
   *
   * @param[in] meshid Mesh id of the degrees of freedom to add the load.
   * @param[in] shapeType Shape function Type (H0 or H1).
   * @param[in] shapeOrder Order of the shape functions to approximate delta u.
   * @param[in] Solution Vector with prescribed solution, length must be a
   * multiple of 3, depending on the length, the load interpolation function is
   * chosen.
   * @param[in] propNumber Proportional load function number for the load.
   * @param[in] direction If local is a function, this vector estimates if the
   * integration is in positive or negative direction of an edge.
   * @param[in] local If true, the loads are treated in the local coordinate
   * system.
   * @param[in] add If true, the loads will be added to the current prescribed
   * values otherwise prescribed values will be overridden.
   */
  void setPrescribedSolution(LoadList &loadlist, indexType meshid,
                             ShapeFunctionTypes shapeType, indexType shapeOrder,
                             Types::VectorX<prec> &Solution,
                             indexType propNumber,
                             Types::VectorX<prec> &direction, bool local,
                             bool add) override;

  /**
   * @brief Creates a std::vector of GenericNodes pointers.
   *
   * @param pointers Pointer to global data.
   * @param[out] nodeVector std::vector containing the nodes of the NodeSet with
   * mesh id "meshId".
   * @param[in] meshId The mesh id of the NodeSet containing the nodes.
   */
  void getNodes(std::vector<GenericNodes *> &nodeVector,
                indexType meshId) override;

  void geometryToParaview(vtkPlotInterface &paraviewAdapter, indexType mainMesh,
                          indexType subMesh);

  void connectEdge(indexType edgeId);
  auto getConnectedEdges() -> std::vector<indexType>;

  void connectFace(indexType faceId);
  auto getConnectedFaces() -> std::vector<indexType>;

  template <typename OStream>
  friend OStream &operator<<(OStream &os, const VertexRuntime &c) {
    fmt::format_to(std::ostream_iterator<char>(os),
                   "Vertex id:  {:>12},  x-coor:  {:>12.4e},  y-coor:  "
                   "{:>12.4e},  z-coor:  {:>12.4e}",
                   c.m_Data_Element->getId(), c.coor(0), c.coor(1), c.coor(2));

    return os;
  }

protected:
  VertexData *m_Vertex_Data_Element;
  Types::Vector3<prec> coor;

private:
  static const HierAMuS::Geometry::GeometryTypes type;
};

} // namespace HierAMuS::Geometry
