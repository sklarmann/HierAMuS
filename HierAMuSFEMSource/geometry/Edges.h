// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "equations/GenericNodes.h"
#include "pointercollection/pointercollection.h"
#include <forwarddeclaration.h>
#include "datatypes.h"
#include <geometry/Base.h>

#include <geometry/GeometryTypes.h>
#include <types/MatrixTypes.h>

#include <vector>

namespace HierAMuS::Geometry {

class Edges : public Base {
public:
  Edges();
  ~Edges() override;
  auto getGroupType() -> const GeometryTypes & override;
  virtual auto getEdgeOrientation(indexType startNode,
                                  indexType endNode) -> prec = 0;

  virtual auto getA1Vector(PointerCollection &pointers,
                   IntegrationPoint &integration_point) -> Types::Vector3<prec> {
    throw std::runtime_error(
        "Function getA1Vector not implemented for the used Edge object!");
    return {};
  };

  virtual auto getJacobian(PointerCollection &pointers, prec xi) -> prec {
    return prec(0);
  };

  auto getJacobian(PointerCollection &pointers,
                           IntegrationPoint &IntegrationPt)
      -> Types::MatrixXX<prec> override {
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
  auto getCoordinates(PointerCollection &pointers, prec xi)
      -> Types::Vector3<prec> override = 0;
  auto getCoordinates(PointerCollection& pointers, IntegrationPoint& IntPoint) -> Types::Vector3<prec> override = 0;

  /**
   * @brief Get the Vertex "number" of object.
   *
   * @param[in] pointers Global data collection.
   * @param[in] number Local Vertex number.
   * @return Geometry::Vertex* Pointer to vertex object.
   */
  virtual auto getVertex(PointerCollection &pointers, indexType number)
      -> Vertex & = 0;

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

  // H1 Shapes
  virtual void setH1Shapes(PointerCollection &pointers, indexType meshId,
                           indexType order, NodeTypes type) {
    Edges::throwError(
        "Error when calling setH1Shapes for Edges, not implemented!");
  };
  virtual void setH1ShapesInternal(PointerCollection &pointers,
                                   indexType meshId, indexType order,
                                   NodeTypes type) {
    Edges::throwError(
        "Error when calling setH1ShapesInternal for Edges, not implemented!");
  };

  virtual void getH1Dofs(PointerCollection &pointers,
                         std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) {
    Edges::throwError(
        "Error when calling getH1Dofs for Edges, not implemented!");
  };

  virtual void getH1DofsInternal(PointerCollection &pointers,
                                 std::vector<DegreeOfFreedom *> &Dofs,
                                 indexType meshID, indexType order) {
    Edges::throwError(
        "Error when calling getH1DofsInternal for Edges, not implemented!");
  };

  virtual auto getH1Nodes(PointerCollection &pointers, indexType meshID,
                          indexType order) -> std::vector<GenericNodes *> {
    Edges::throwError(
        "Error when calling getH1Nodes for Edges, not implemented!");
    return {};
  };

  virtual auto getH1NodesInternal(PointerCollection &pointers, indexType meshID,
                                  indexType order)
      -> std::vector<GenericNodes *> {
    Edges::throwError(
        "Error when calling getH1Nodes for Edges, not implemented!");
    return {};
  };

  virtual void getH1Shapes(PointerCollection &pointers, indexType order,
                           Types::VectorX<prec> &shape,
                           Types::VectorX<prec> &shapeDerivative, prec xsi) {
    Edges::throwError(
        "Error when calling getH1Shapes for Edges, not implemented!");
  };

  virtual void getH1ShapesInternal(PointerCollection &pointers, indexType order,
                                   Types::VectorX<prec> &shape,
                                   Types::VectorX<prec> &shapeDerivative,
                                   prec xsi) {
    Edges::throwError(
        "Error when calling getH1ShapesInternal for Edges, not implemented!");
  };

  virtual auto getH1Shapes(PointerCollection &pointers, indexType order,
                           IntegrationPoint &integration_point) -> H1Shapes {
    Edges::throwError(
        "Error when calling getH1Shapes for Edges, not implemented!");
    return {};
  };

  virtual auto getH1ShapesInternal(PointerCollection &pointers, indexType order,
                                   IntegrationPoint &integration_point)
      -> H1Shapes {
    Edges::throwError(
        "Error when calling getH1Shapes for Edges, not implemented!");
    return {};
  };

  // HDivShapes
  void setHDivShapes(PointerCollection &pointers, indexType meshid,
                     indexType order, NodeTypes type) {
    Edges::throwError(
        "Error when calling setHDivShapes for Faces, not implemented!");
  };

  void getHDivDofs(PointerCollection &pointers,
                   std::vector<DegreeOfFreedom *> &Dofs,
                   indexType meshID, indexType order,
                   NodeTypes type) {
    Edges::throwError(
        "Error when calling getHDivDofs for Faces, not implemented!");
  };
  void getHDivShapes(PointerCollection &pointers, indexType order,
                     Types::Matrix22<prec> &jacobi, Types::VectorX<prec> &shape,
                     Types::Matrix2X<prec> &dshape, prec xi,
                     prec eta) {
    HierAMuS::Geometry::Edges::throwError(
        "Error when calling getHDivShapes for Faces, not implemented!");
  };

  virtual auto getDirectionVector(PointerCollection &pointers) -> Types::Vector3<prec> {
    return Types::Vector3<prec>(0);
  };

  void setAllNodeBoundaryConditionMeshId(PointerCollection &pointers,
                                                 indexType meshId,
                                                 indexType dof) override;

  /**
   * @brief Changes orientation of the face (reverts the face normal vector).
   */
  virtual void flip() = 0;

private:
  static const HierAMuS::Geometry::GeometryTypes type;
  static void throwError(const std::string &msg) { throw std::runtime_error(msg); };
};

} // namespace HierAMuS
