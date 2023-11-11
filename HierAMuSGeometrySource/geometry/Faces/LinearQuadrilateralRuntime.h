// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <array>

#include "geometry/VertexRuntime.h"

#include "geometry/Faces/FacesH1Interface.h"
#include "geometry/Faces/FacesHDivInterface.h"
#include "geometry/Faces/FacesRuntimeDataInterface.h"

#include <types/MatrixTypes.h>

namespace HierAMuS::Geometry {
class LinearQuadrilateralData;
/**
 * @brief Linear Quadrilateral Geometry Element. Represents a 4-vertex face
 * element.
 *
 */
class LinearQuadrilateralRuntime
    : public FacesRuntimeDataInterface<4, 4, LinearQuadrilateralData,
                                       LinearQuadrilateralRuntime>,
      public FacesH1Interface,
      public FacesHDivInterface {
public:
  LinearQuadrilateralRuntime(GeometryData &geoData,
                             LinearQuadrilateralData &base_element);
  ~LinearQuadrilateralRuntime() override;

  auto getH1Face() -> FacesH1Interface * override { return this; };
  auto getHDivFace() -> FacesHDivInterface * override { return this; };

  /**
   * @brief Get the Type object
   *
   * @return const GeometryTypes&, type of the element, here
   * GeometryTypes::LinearQuadrilateral
   */
  auto getType() -> const GeometryTypes & override;

  /**
   * @brief Prints the object to the console and to file.
   *
   * @param pointers, object containing the pointers to global data.
   */
  void print(spdlog::logger &Log) override;

  /**
   * @brief Get the tangent vector G1 of the element at xi, eta, G1 is X,xi.
   *
   * @param pointers object containing the pointers to global data.
   * @param integrationPoint integration point.
   * @return Types::Vector3<prec> tangent vector G1 of the element at xi, eta.
   */
  auto getTangent_G1(IntegrationPoint &integrationPoint)
      -> Types::Vector3<prec> override;
  /**
   * @brief Get the tangent vector G2 of the element at xi, eta, G2 is X,eta.
   *
   * @param pointers object containing the pointers to global data.
   * @param integrationPoint integration point.
   * @return Types::Vector3<prec> tangent vector G2 of the element at xi, eta.
   */
  auto getTangent_G2(IntegrationPoint &integrationPoint)
      -> Types::Vector3<prec> override;

  /**
   * @brief Get the Face Normal vector
   *
   * @param pointers[in], pointers to global data.
   * @return Types::Vector3<prec>, Normal vector of the face
   */
  auto getFaceNormal()
      -> Types::Vector3<prec> override;

  /**
   * @brief Get the Orientation object
   *
   * @param vertex1 Vertex 1 of comparison edge
   * @param vertex2 Vertex 2 of comparison edge
   * @return faceorientation
   */
  auto getOrientation(indexType vertex1,
                      indexType vertex2) -> faceorientation override;

  /**
   * @brief Modifies the current integration point and evaluates a shape factor
   * based on the orientation.
   *
   * @param[inout] IP Current integration point.
   * @param[out] shapeFactor factor to scale the shape function.
   * @param[out] orientation orientation of the face.
   */
  void modifyIntegrationpoint(IntegrationPoint &IP, prec &shapeFactor,
                              faceorientation orientation) override;

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
  auto getH1Dofs(indexType meshID, indexType order)
      -> std::vector<DegreeOfFreedom *> override;
  void getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                 indexType order) override;
  void getH1DofsInternal(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) override;
  auto getH1Nodes(indexType meshID, indexType order)
      -> std::vector<GenericNodes *> override;

  auto getH1NodesInternal(indexType meshID, indexType order)
      -> std::vector<GenericNodes *> override;

  auto getH1NodesList(indexType meshID,
                      indexType order) -> MeshIdNodeList override;

  auto getHDivNodes(indexType meshID,
                    indexType order) -> std::vector<GenericNodes *> override;

  auto getHDivNodesInternal(indexType meshID,
                            indexType order)
      -> std::vector<GenericNodes *> override;

  /**
   * @brief New version of getH1Shapes.
   * Computes the Edge and Vertex H1 shapes for the quadrilateral face.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param order[in], order of the H1 shape functions.
   * @param IntegrationPt[in], integration point.
   * @return H1Shapes, H1 shape functions.
   */
  auto getH1Shapes(indexType order, IntegrationPoint &IntegrationPt)
      -> H1Shapes override;
  /**
   * @brief New version of getH1ShapesInternal.
   * Computes the internal H1 shape functions known as bubble functions.
   *
   * @param[in] pointers object containing the pointers to global data.
   * @param[in] order order of the shape functions.
   * @param[in] IntegrationPt Current integration point.
   * @param[in] orientation Orientation of the face (used for bubble functions).
   * @return H1Shapes Object containing the H1 shape functions.
   */
  auto getH1ShapesInternal(indexType order, IntegrationPoint &IntegrationPt,
                           faceorientation orientation = faceorientation::p_1)
      -> H1Shapes override;

  // HDivShapes
  void setHDivShapes(indexType meshid,
                     indexType order, NodeTypes type) override;

  void getHDivDofs(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                   indexType order, NodeTypes type) override;
  void getHDivShapes(indexType order,
                     Types::Matrix2X<prec> &shape, Types::VectorX<prec> &dshape,
                     prec xi, prec eta) override;

  auto getHDivShapes(indexType order,
                     IntegrationPoint &IntegrationPt) -> HDivShapes override;

  // Curl Shapes

  // Special Plate Shapes
  void setSpecialPlateShapes(indexType meshid,
                             indexType order, NodeTypes type) override;
  auto getSpecialPlateDofs(indexType meshID,
                           indexType order, NodeTypes type)
      -> std::vector<DegreeOfFreedom *> override;
  auto getSpecialPlateShapes(IntegrationPoint &intPoint, indexType order)
      -> SpecialPlateShapes override;

  /**
   * @brief Sets boundary condition at degree of freedom "dof" on all nodes in
   * NodeSet with mesh id "meshid".
   *
   * @param pointers Pointers to global data.
   * @param meshId Mesh id of the NodeSet.
   * @param dof Degree of freedom number to set the boundary condition, zero
   * based indexing.
   */
  void setAllNodeBoundaryConditionMeshId(indexType meshId,
                                         indexType dof) override;

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
  void setBoundaryCondition(indexType meshId, indexType order,
                            ShapeFunctionTypes shapeType,
                            Types::Vector3<indexType> &dofs, bool set) override;

  // Paraview part
  /**
   * @brief Transfers the geometry of the element to the Paraview Adapter
   * (vtkPlotInterface).
   *
   * @param paraviewAdapter[in], ParaviewAdapter object.
   * @param mainMesh [in], main mesh, normally the value 0.
   * @param subMesh[in], sub mesh normally the material number.
   */
  void geometryToParaview(vtkPlotInterface &paraviewAdapter, indexType mainMesh,
                          indexType subMesh) override;
  /**
   * @brief Compute weights for the vertices on the Paraview mesh. Required when
   * stress (or other values only available a the integration points) projection
   * is done.
   *
   * @param paraviewAdapter[in], ParaviewAdapter object.
   * @param mainMesh[in], main mesh, normally the value 0.
   * @param subMesh[in], sub mesh normally the material number.
   */
  void computeWeightsParaview(vtkPlotInterface &paraviewAdapter,
                              indexType mainMesh, indexType subMesh) override;
  /**
   * @brief Transfers a solution field to the Paraview Adapter
   * (vtkPlotInterface).
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param paraviewAdapter[in], ParaviewAdapter object.
   * @param mainMesh[in], main mesh, normally the value 0.
   * @param subMesh[in], sub mesh normally the material number.
   * @param meshId[in], mesh id of the solution field which is transferred.
   * @param order[in], order of the approximation of the solution field.
   * @param name[in], name of the solution field, predefined names are in
   * paraviewNames.
   */
  void H1SolutionToParaview(vtkPlotInterface &paraviewAdapter,
                            indexType mainMesh, indexType subMesh,
                            indexType order, Types::VectorX<prec> &solution,
                            std::string &name) override;
  /**
   * @brief Transfers H1 data of given order to the vertices of the Paraview
   * Adapter. Currently transfers only the linear approximation of the data.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param paraviewAdapter[in], ParaviewAdapter object.
   * @param mainMesh[in], main mesh, normally the value 0.
   * @param subMesh[in], sub mesh normally the material number.
   * @param Data[in], data to be transferred, Vector with the size
   * numberComponents * numberOfNodes to represent the approximation order.
   * @param numberComponents[in], number of components of the data.
   * @param order[in], order of the approximation of the data.
   * @param name[in], name of the data, predefined names are in paraviewNames.
   */
  void H1DataToParaview(vtkPlotInterface &paraviewAdapter, indexType mainMesh,
                        indexType subMesh, Types::VectorX<prec> &Data,
                        indexType numberComponents, indexType order,
                        std::string &name) override;
  /**
   * @brief Projects the the data at the current Integration point to the
   * vertices in the Paraview mesh.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param paraviewAdapter[in], ParaviewAdapter object.
   * @param mainMesh[in], main mesh, normally the value 0.
   * @param subMesh[in], submesh normally the material number.
   * @param order[in], order of the approximation of the data.
   * @param IntegrationPt[in], Integration point to be projected.
   * @param data[in], data to be projected.
   * @param numberComponents[in], number of components of the data.
   * @param name[in], name of the data, predefined names are in paraviewNames.
   */
  void projectDataToParaviewVertices(
      vtkPlotInterface &paraviewAdapter,
      indexType mainMesh, indexType subMesh, indexType order,
      IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
      indexType numberComponents, std::string name) override;

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

  void set_geometry_pointers(GeometryData &geoData) override;

private:
  static const GeometryTypes type;

protected:
};

} // namespace HierAMuS::Geometry
