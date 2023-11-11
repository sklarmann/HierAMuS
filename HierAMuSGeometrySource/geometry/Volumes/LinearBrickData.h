// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once


#include <geometry/GeometryTypes.h>
#include <geometry/Volumes/VolumesDataInterface.h>

#include <types/MatrixTypes.h>

#include <vector>

namespace HierAMuS::Geometry {

/**
 * @brief Linear Brick geometry element. Represents a 8-vertex, 12-edges,
 * 6-faces volume element. The local numbering of the geometric subelements
 * (vertices, edges, faces) is counterclockwise from bottom to top.
 * @image html geometryelements/linearbricknumbering.svg
 *
 * The vertex numbers in the given picture are black, the edge numbers are blue,
 * and the face numbers are red. The numbers given, indicate the order in which
 * the vertex, edge, and face numbers need to be passed to the geometric
 * element.
 *
 */
class LinearBrickData : public VolumesDataInterface<8,12,6,LinearBrickData> {
public:
  LinearBrickData();
  ~LinearBrickData() override;
  auto getType() -> const GeometryTypes & override;

  auto getRuntimeObject(GeometryData &geoData)
      -> std::shared_ptr<VolumesRuntime> override;

  /**
   * @brief Get the Integration Points for the brick element.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @return IntegrationPoints
   */
  auto getIntegrationPoints(indexType elementId) -> IntegrationPoints override;

  
  // H1 Shapes
  /**
   * @brief Assigns the nodes to the geometry object for H1 shape functions
   * with given order associated with the specific meshId. A single element can
   * have multiple H1 nodes with different mesh ids (e.g. for different solutions
   * like displacements and rotations).
   *
   * @param pointers, object containing the pointers to global data.
   * @param meshId, mesh id of the nodes for the H1 shape functions.
   * @param order, order of the H1 shape functions.
   * @param type, type of the degrees of freedom, currently unused.
   */
  void setH1Shapes(indexType meshId,
                   indexType order, NodeTypes type) override;
  void setH1ShapesInternal(indexType meshId,
                           indexType order, NodeTypes type) override;

  
  
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
  void checkUpdateElement(EquationHandler &eqHandler,
                          GeometryData &geoData) override;


  void setAllNodeBoundaryConditionMeshId(indexType meshId, indexType dof) override;

  
  
  static std::string getName() { return "Linear Brick"; };

private:
  
  static const GeometryTypes type;
};

} // namespace HierAMuS