// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "geometry/Edges/EdgesDataInterface.h"

#include "geometry/GeometryTypes.h"
#include <types/MatrixTypes.h>

#include <array>
#include <vector>

namespace HierAMuS::Geometry {

class LinearEdgeData : public EdgesDataInterface<2, LinearEdgeData> {
public:
  LinearEdgeData();
  ~LinearEdgeData() override;

  auto getRuntimeObject(GeometryData &geoData)
      -> std::shared_ptr<EdgesRuntime> override;

  auto getType() -> const GeometryTypes & override;

  void getNodes(std::vector<GenericNodes *> &nodeVector,
                indexType meshId) override;

  auto getH1NodesList(indexType meshID,
                      indexType order) -> MeshIdNodeList override;
  /**
   * @brief Creates a std::vector of GenericNodes pointers of Internal nodes.
   *
   * @param pointers Pointer to global data.
   * @param[out] nodeVector std::vector containing the nodes of the NodeSet with
   * mesh id "meshId".
   * @param[in] meshId The mesh id of the NodeSet containing the nodes.
   */
  void getNodesInternal(std::vector<GenericNodes *> &nodeVector,
                        indexType meshId) override;

  auto getEdgeOrientation(indexType startVertex, indexType endVertex)
      -> prec override;
  /**
   * @brief Set Boundary conditions on the element
   *
   * @param[in] pointers Pointers to global data.
   * @param[in] meshId Mesh id of the nodes.
   * @param[in] order Order of the shape function.
   * @param[in] shapeType Type of the shape function to restrict.
   * @param[in] dofs Degree of freedom list to set boundary conditions.
   * @param[in] set If true, boundary conditions are overridden, otherwise only
   * 1 is considered.
   */
  void setBoundaryCondition(indexType meshId,
                            indexType order, ShapeFunctionTypes shapeType,
                            Types::Vector3<indexType> &dofs, bool set) override;

  /**
   * @brief Computes the physical coordinates at local coordinate xi.
   *
   * @param pointers Global data collection.
   * @param xi Local coordinate xi.
   * @return Returns a vector with the x, y, z coordinates.
   */
  auto getCoordinates(prec xi)
      -> Types::Vector3<prec> override;
  auto getCoordinates(IntegrationPoint &IntPoint)
      -> Types::Vector3<prec> override;

  /**
   * @brief Get the Integration Points object. Returns the IntegrationPoints
   * object, already set to the correct type. Only the integration order needs
   * to be specified.
   *
   * @param[in] pointers, object containing the pointers to global data.
   * @return IntegrationPoints, IntegrationPoints object set to the correct
   * type.
   */
  auto getIntegrationPoints(indexType elementId) -> IntegrationPoints override;

  // H1Shapes
  void setH1Shapes(indexType meshId,
                   indexType order, NodeTypes type) override;
  void setH1ShapesInternal(indexType meshId,
                           indexType order, NodeTypes type) override;

  void getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                 indexType order) override;
  void getH1DofsInternal(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) override;

  auto getH1Nodes(indexType meshID,
                  indexType order) -> std::vector<GenericNodes *> override;

  auto getH1NodesInternal(indexType meshID,
                          indexType order)
      -> std::vector<GenericNodes *> override;
  /**
   * @brief Get Shape functions and derivatives with respect to parameter space.
   *
   * @tparam prec
   * @tparam indexType
   * @param[in] pointers Collection of global data.
   * @param[in] order Order of the shape functions.
   * @param[out] shape Shape functions.
   * @param[out] shapeDerivative Derivatives of shape functions with respect to
   * parameter space.
   * @param[out] xsi Position in parameter space, ranging from -1 to 1.
   */
  void getH1Shapes(indexType order,
                   Types::VectorX<prec> &shape,
                   Types::VectorX<prec> &shapeDerivative, prec xsi) override;
  /**
   * @brief Routine for the internal shape functions of the edge.
   *
   * @tparam prec
   * @tparam indexType
   * @param[in] pointers Collection of global data.
   * @param[in] order Order of the shape functions.
   * @param[out] shape Shape functions.
   * @param[out] shapeDerivative Derivatives of shape functions with respect to
   * parameter space.
   * @param[in] xsi Position in parameter space, ranging from -1 to 1.
   */
  void getH1ShapesInternal(indexType order,
                           Types::VectorX<prec> &shape,
                           Types::VectorX<prec> &shapeDerivative,
                           prec xsi) override;

  /**
   * @brief Routine for the shape functions of the edge.
   *
   * @tparam prec
   * @tparam indexType
   * @param[in] pointers Collection of global data.
   * @param[in] order Order of the shape functions.
   * @param[in] integration_point Current integration point.
   * @return H1Shapes, Struct containing the shape functions and derivatives in
   * parameter space.
   */
  auto getH1Shapes(indexType order,
                   IntegrationPoint &integration_point) -> H1Shapes override;
  /**
   * @brief Routine for the internal shape functions of the edge.
   *
   * @tparam prec
   * @tparam indexType
   * @param[in] pointers Collection of global data.
   * @param[in] order Order of the shape functions.
   * @param[in] integration_point Current integration point.
   * @return H1Shapes, Struct containing the shape functions and derivatives in
   * parameter space.
   */
  auto getH1ShapesInternal(indexType order,
                           IntegrationPoint &integration_point)
      -> H1Shapes override;

  auto getJacobian(prec xi) -> prec override;
  auto getJacobian(IntegrationPoint &IntegrationPt)
      -> prec override;

  /**
   * @brief Sets load using a Geometric Object.
   *
   * @param[in] pointers Pointers to global data.
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
   * @param[in] pointers Pointers to global data.
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

  

  auto getDirectionVector()
      -> Types::Vector3<prec> override;

  auto getA1Vector(IntegrationPoint &integration_point)
      -> Types::Vector3<prec> override;

  void setAllNodeBoundaryConditionMeshId(indexType meshId,
                                         indexType dof) override;

  /**
   * @brief Changes orientation of the face (reverts the face normal vector).
   */
  void flip() override;

  static auto getName() -> std::string { return "Linear Edge"; };

private:
  static const GeometryTypes type;
};

} // namespace HierAMuS::Geometry
