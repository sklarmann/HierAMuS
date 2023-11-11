// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Nodetypes.h"
#include <geometry/GeometryTypes.h>
#include <shapefunctions/IntegrationsPoints/IntegrationPoints.h>
#include <stdexcept>
#include <types/MatrixTypes.h>
#include <vector>

#include "datatypes.h"
#include "NodeSetNodeList.h"
#include "NodeSetList.h"
#include "NodeSetManager.h"

#include <tuple>

#include <iostream>

namespace HierAMuS {
class NodeSet;
class vtkPlotInterface;
class GenericNodes;
class LoadList;

namespace Geometry {


class GeometryData;

/**
 * @brief Base geometry class.
 * All other geometry classes are derived from this class.
 * This class handles the connection to the degrees of freedom of a geometric
 * object.
 *
 */
class GeometryBaseData {

public:
  GeometryBaseData();
  virtual ~GeometryBaseData();

  /**
   * @brief Set the id of the geometric object for identification.
   *
   * @param[in] id Id of the object.
   */
  void set_id(indexType id) { this->id = id; };
  /**
   * @brief Get the id of the current object.
   *
   * @return Returns the if of the object of type indexType.
   */
  auto getId() -> indexType { return this->id; };

  /**
   * @brief Set Boundary conditions on the element
   *
   * @param[in] pointers Pointers to global data.
   * @param[in] meshId Mesh id of the nodes.
   * @param[in] order Order of the shape function.
   * @param[in] shapeType Type of the shape function to restrict.
   * @param[in] dofs Degree of freedom list to set boundary conditions.
   * @param[in] set If true, boundary conditions are overridden, otherwise only 1
   * is considered.
   */
  virtual void setBoundaryCondition(indexType meshId,
                       indexType order, ShapeFunctionTypes shapeType,
                       Types::Vector3<indexType> &dofs,
                       bool set){}; // move into specialized classes

  /**
   * @brief Set a NodeSet for the geometric element containing the nodes.
   *
   * @param pointers Global data collection.
   * @param[in] meshID The mesh id of the nodes.
   * @param[in] numberOfNodes Number of nodes the NodeSet contains.
   * @param[in] type The type of the nodes.
   */
  void setNodeSet(indexType meshID,
                  indexType numberOfNodes, NodeTypes type);


  /**
   * @brief Creates a std::vector of GenericNodes pointers.
   *
   * @param pointers Pointer to global data.
   * @param[out] nodeVector std::vector containing the nodes of the NodeSet with
   * mesh id "meshId".
   * @param[in] meshId The mesh id of the NodeSet containing the nodes.
   */
  auto getNodesOfSet(indexType meshId) -> std::vector<GenericNodes *>;

  auto getNodeSetNodeListMeshId(indexType meshId) -> NodeSetNodeList;

  auto getNodeSetList() -> NodeSetList;

  /**
   * @brief Creates a std::vector of GenericNodes pointers.
   *
   * @param pointers Pointer to global data.
   * @param[out] nodeVector std::vector containing the nodes of the NodeSet with
   * mesh id "meshId".
   * @param[in] meshId The mesh id of the NodeSet containing the nodes.
   */
  virtual void getNodes(std::vector<GenericNodes *> &nodeVector,
                        indexType meshId){}; // move into specialized classes
  /**
   * @brief Creates a std::vector of GenericNodes pointers of Internal nodes.
   *
   * @param pointers Pointer to global data.
   * @param[out] nodeVector std::vector containing the nodes of the NodeSet with
   * mesh id "meshId".
   * @param[in] meshId The mesh id of the NodeSet containing the nodes.
   */
  virtual void
  getNodesInternal(std::vector<GenericNodes *> &nodeVector,
                   indexType meshId); // move into specialized classes

  /**
   * @brief Sets boundary condition at degree of freedom "dof" on all nodes in
   * NodeSet with mesh id "meshid".
   *
   * @param pointers Pointers to global data.
   * @param meshId Mesh id of the NodeSet.
   * @param dof Degree of freedom number to set the boundary condition, zero
   * based indexing.
   */
  virtual void setAllNodeBoundaryConditionMeshId(indexType meshId,
                                                 indexType dof);
  void getAllEquationsIds(std::vector<DegreeOfFreedom *> &Dofs);
  void getNodeEquationIds(std::vector<DegreeOfFreedom *> &Dofs,
                          indexType meshId, indexType nodeNumber);

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
  virtual void setLoad(LoadList &loadlist, indexType meshid,
                       ShapeFunctionTypes shapeType, indexType shapeOrder,
                       Types::VectorX<prec> &Loads, indexType propNumber,
                       Types::VectorX<prec> &direction, bool local,
                       bool add){}; // move into specialized classes

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
  virtual void
  setPrescribedSolution(LoadList &loadlist, indexType meshid,
                        ShapeFunctionTypes shapeType, indexType shapeOrder,
                        Types::VectorX<prec> &Solution, indexType propNumber,
                        Types::VectorX<prec> &direction, bool local,
                        bool add){}; // move into specialized classes

  /**
   * @brief Returns the type of the geometric object.
   *
   * @return Type of the object.
   */
  virtual auto getType() -> const GeometryTypes &;
  /**
   * @brief Returns the group type of the geometric element.
   * The group types are: Edges, Faces, Volumes.
   *
   * @return GeometryTypes as the type
   */
  virtual auto getGroupType() -> const GeometryTypes &;

  // Output Information
  virtual void print(spdlog::logger &Logger) = 0;
  void printEqInfo(spdlog::logger &Logger);

  // Getters and setters



  auto getNodeListMap() -> std::map<indexType, NodeSetNodeList>;
  



  virtual auto getIntegrationPoints(indexType elementId) -> IntegrationPoints; // move into specialized classes

  void setNodeSetManager(NodeSetManager &other) { m_NodeSetManager = other; };
  void setNodeSetManager(NodeSetManager &&other) { m_NodeSetManager = other; };

protected:
  indexType id;

  NodeSetManager m_NodeSetManager;

private:
  static const HierAMuS::Geometry::GeometryTypes type;
};

} // namespace Geometry
} // namespace HierAMuS
