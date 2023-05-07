// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "pointercollection/pointercollection.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <array>
#include <forwarddeclaration.h>

#include <geometry/Faces.h>
#include <types/MatrixTypes.h>

namespace HierAMuS::Geometry {

/**
 * @brief Linear Quadrilateral Geometry Element. Represents a 4-vertex face
 * element.
 *
 */
class LinearQuadrilateral : public Faces {
public:
  LinearQuadrilateral();
  ~LinearQuadrilateral() override;
  /**
   * @brief Get the Type object
   *
   * @return const GeometryTypes&, type of the element, here
   * GeometryTypes::LinearQuadrilateral
   */
  auto getType() -> const GeometryTypes & override;
  /**
   * @brief Get the Number of Vertices object
   *
   * @return indexType, number of vertices of the element, here 4
   */
  auto getNumberOfVerts() -> indexType override { return 4; };
  /**
   * @brief Returns the number of edges of the element.
   * @return Number of edges as an indexType
   */
  auto getNumberOfEdges() -> indexType override { return 4; };
  /**
   * @brief Set the Vertex numbers of the object.
   *
   * @param[in] vertsIn, vector of 4 integers of type indexType refering to the
   * associated vertex numbers.
   */
  void setVerts(GeometryData &geoData, std::vector<indexType> &vertsIn) override;
  /**
   * @brief Set the Edges of the object. In total the element has 4 edges.
   *
   * @param[in] edgesIn, vector of 4 integers of type indexType refering to the
   * associated edge numbers.
   */
  void setEdges(const std::vector<indexType> &edgesIn) override;
  /**
   * @brief Get the global vertex numbers of the object.
   *
   * @param[out] vertsOut, vector of 4 integers of type indexType with global
   * vertex numbers of the element.
   */
  void getVerts(std::vector<indexType> &vertsOut) override;
  /**
   * @brief Get the Verts object
   *
   * @param[in] pointers, object containtig the pointers to global data.
   * @param[out] vertsOut, vector of Geometry::Base pointers to the vertices of
   * the geometry element.
   */
  void getVerts(PointerCollection &pointers,
                std::vector<Base *> &vertsOut) override;
  /**
   * @brief Get a pointer to the local vertex local_number of the quadrilateral
   * face.
   *
   * @param pointers object containtig the pointers to global data.
   * @param local_number local number of the vertex.
   * @return Geometry::Vertex* pointer to the vertex.
   */
  auto getVertex(PointerCollection &pointers, indexType local_number)
      -> Geometry::Vertex * override;
  /**
   * @brief Get a pointer to the local edge local_number of the quadrilateral
   * face.
   *
   * @param pointers object containtig the pointers to global data.
   * @param local_number local number of the edge.
   * @return Geometry::Vertex* pointer to the edge.
   */
  auto getEdge(PointerCollection &pointers, indexType local_number)
      -> Geometry::Edges * override;
  /**
   * @brief Get the global edge numbers of the object.
   *
   * @param[out] edgesOut, vector of 4 integers of type indexType with global
   * edge numbers of the element.
   */
  void getEdges(std::vector<indexType> &edgesOut) override;
  /**
   * @brief Get the Edges object
   *
   * @param pointers, object containtig the pointers to global data.
   * @param edgesOut, vector of Geometry::Base pointers to the edges of the
   * geometry element.
   */
  void getEdges(PointerCollection &pointers,
                std::vector<Base *> &edgesOut) override;

  /** @brief Returns the global edge number based on the local edge number of the face.
   *
   *  detailed description
   *
   *  @param[in] local_number local edge number of the face
   *  @return indexType global edge number
   */
  auto getEdge(indexType local_number) -> indexType override;

  /** @brief Check if the face has the three vertices.
    *
    * @param[in] v1 Global number of a corner vertex.
    * @param[in] v2 Global number of a corner vertex.
    * @param[in] v3 Global number of a corner vertex.
    * @return true If the face has the three vertices.
    * @return false If the face does not have the three vertices.
    */
  auto hasVertices(indexType v1, indexType v2, indexType v3) -> bool override;
  /**
   * @brief Prints the object to the console and to file.
   *
   * @param pointers, object containtig the pointers to global data.
   */
  void print(PointerCollection &pointers) override;
  /**
   * @brief Get the Coordinates x,y,z at xi, eta of the current element.
   *
   * @param[in] pointers, object containtig the pointers to global data.
   * @param[in] xi, local xi coordinate.
   * @param[in] eta, local eta coordinate.
   * @return Types::Vector3<prec>, global coordinates x,y,z of the element at
   * xi, eta.
   */
  auto getCoordinates(PointerCollection &pointers, prec xi, prec eta)
      -> Types::Vector3<prec> override;

  auto getCoordinates(PointerCollection &pointers,
                      IntegrationPoint &integrationPoint)
      -> Types::Vector3<prec> override;

  /**
   * @brief Get the tangent vector G1 of the element at xi, eta, G1 is X,xi.
   *
   * @param pointers object containtig the pointers to global data.
   * @param integrationPoint integration point.
   * @return Types::Vector3<prec> tangent vector G1 of the element at xi, eta.
   */
  auto getTangent_G1(PointerCollection &pointers,
                     IntegrationPoint &integrationPoint)
      -> Types::Vector3<prec> override;
  /**
   * @brief Get the tangent vector G2 of the element at xi, eta, G2 is X,eta.
   *
   * @param pointers object containtig the pointers to global data.
   * @param integrationPoint integration point.
   * @return Types::Vector3<prec> tangent vector G2 of the element at xi, eta.
   */
  auto getTangent_G2(PointerCollection &pointers,
                     IntegrationPoint &integrationPoint)
      -> Types::Vector3<prec> override;

  /**
   * @brief Get the Face Normal vector
   *
   * @param pointers[in], pointers to global data.
   * @return Types::Vector3<prec>, Normal vector of the face
   */
  auto getFaceNormal(PointerCollection &pointers)
      -> Types::Vector3<prec> override;

  /**
   * @brief Get the Orientation object
   *
   * @param pointers PointerCollection object containtig the pointers to global
   * data.
   * @param vertex1 Vertex 1 of comparison edge
   * @param vertex2 Vertex 2 of comparison edge
   * @return faceorientation
   */
  auto getOrientation(PointerCollection &pointers, indexType vertex1,
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
   * @param[in] pointers, object containtig the pointers to global data.
   * @return IntegrationPoints, IntegrationPoints object set to the correct
   * type.
   */
  auto getIntegrationPoints(PointerCollection &pointers, indexType elementId)
  -> IntegrationPoints override;

  // Jacobian matrix
  /**
   * @brief Get the Jacobian matrix at the current IntegrationPoint. New version
   *
   * @param pointers, object containtig the pointers to global data.
   * @param IntegrationPt, current IntegrationPoint.
   * @return Types::MatrixXX<prec>, Jacobian matrix at the current
   * IntegrationPoint.
   */
  auto getJacobian(PointerCollection &pointers, IntegrationPoint &IntegrationPt)
      -> Types::MatrixXX<prec> override;
  // H1Shapes
  /**
   * @brief Assignes the nodes to the geometry object for H1 shape functions
   * with given order associated with the specific meshId. A single element can
   * have multiple H1 nodes with different meshids (e.g. for different solutions
   * like displacements and rotations).
   *
   * @param pointers, object containtig the pointers to global data.
   * @param meshId, mesh id of the nodes for the H1 shape functions.
   * @param order, order of the H1 shape functions.
   * @param type, type of the degrees of freedom, currently unused.
   */
  void setH1Shapes(PointerCollection &pointers, indexType meshId,
                   indexType order, NodeTypes type) override;
  void setH1ShapesInternal(PointerCollection &pointers, indexType meshId,
                           indexType order, NodeTypes type) override;
  auto getH1Dofs(PointerCollection &pointers, indexType meshID, indexType order)
      -> std::vector<DegreeOfFreedom *> override;
  void getH1Dofs(PointerCollection &pointers,
                 std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                 indexType order) override;
  void getH1DofsInternal(PointerCollection &pointers,
                         std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) override;
  auto getH1Nodes(PointerCollection &pointers, indexType meshID,
                  indexType order) -> std::vector<GenericNodes *> override;

  auto getH1NodesInternal(PointerCollection &pointers, indexType meshID,
                          indexType order)
      -> std::vector<GenericNodes *> override;

  auto getHDivNodes(PointerCollection &pointers, indexType meshID,
                    indexType order) -> std::vector<GenericNodes *> override;

  auto getHDivNodesInternal(PointerCollection &pointers, indexType meshID,
                            indexType order)
      -> std::vector<GenericNodes *> override;

  void getH1Shapes(PointerCollection &pointers, indexType order,
                   Types::VectorX<prec> &shape,
                   Types::Matrix2X<prec> &shapeDerivative, prec xsi,
                   prec eta) override;

  void getH1ShapesInternal(PointerCollection &pointers, indexType order,
                           Types::VectorX<prec> &shape,
                           Types::Matrix2X<prec> &shapeDerivative, prec xsi,
                           prec eta) override;

  /**
   * @brief New version of getH1Shapes.
   * Computes the Edge and Vertex H1 shapes for the quadrilateral face.
   *
   * @param pointers[in], object containtig the pointers to global data.
   * @param order[in], order of the H1 shape functions.
   * @param IntegrationPt[in], integration point.
   * @return H1Shapes, H1 shape functions.
   */
  auto getH1Shapes(PointerCollection &pointers, indexType order,
                   IntegrationPoint &IntegrationPt) -> H1Shapes override;
  /**
   * @brief New version of getH1ShapesInternal.
   * Computes the internal H1 shape functions known as bubble functions.
   *
   * @param[in] pointers object containtig the pointers to global data.
   * @param[in] order order of the shape functions.
   * @param[in] IntegrationPt Current integration point.
   * @param[in] orientation Orientation of the face (used for bubble functions).
   * @return H1Shapes Object containing the H1 shape functions.
   */
  auto getH1ShapesInternal(PointerCollection &pointers, indexType order,
                           IntegrationPoint &IntegrationPt,
                           faceorientation orientation = faceorientation::p_1)
      -> H1Shapes override;

  // HDivShapes
  void setHDivShapes(PointerCollection &pointers, indexType meshid,
                     indexType order, NodeTypes type) override;

  void getHDivDofs(PointerCollection &pointers,
                   std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                   indexType order, NodeTypes type) override;
  void getHDivShapes(PointerCollection &pointers, indexType order,
                     Types::Matrix2X<prec> &shape, Types::VectorX<prec> &dshape,
                     prec xi, prec eta) override;

  auto getHDivShapes(PointerCollection &pointers, indexType order,
                     IntegrationPoint &IntegrationPt) -> HDivShapes override;

  // Curl Shapes

  // Special Plate Shapes
  void setSpecialPlateShapes(PointerCollection &pointers, indexType meshid,
                             indexType order, NodeTypes type) override;
  auto getSpecialPlateDofs(PointerCollection &pointers, indexType meshID,
                           indexType order, NodeTypes type)
      -> std::vector<DegreeOfFreedom *> override;
  auto getSpecialPlateShapes(PointerCollection &pointers,
                             IntegrationPoint &intPoint, indexType order)
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
  void setAllNodeBoundaryConditionMeshId(PointerCollection &pointers,
                                         indexType meshId,
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
  void setLoad(PointerCollection &pointers, indexType meshid,
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
  void setPrescribedSolution(PointerCollection &pointers, indexType meshid,
                             ShapeFunctionTypes shapeType, indexType shapeOrder,
                             Types::VectorX<prec> &Solution,
                             indexType propNumber,
                             Types::VectorX<prec> &direction, bool local,
                             bool add) override;
  /**
   * @brief Set Boundary conditions on the element
   *
   * @param[in] pointers Pointers to global data.
   * @param[in] meshId Mesh id of the nodes.
   * @param[in] order Order of the shape function.
   * @param[in] shapeType Type of the shape function to restrict.
   * @param[in] dofs Degree of freedom list to set boundary conditions.
   * @param[in] set If true, boudary contidions are overridden, otherwise only 1
   * is considered.
   */
  void setBoundaryCondition(PointerCollection &pointers, indexType meshId,
                            indexType order, ShapeFunctionTypes shapeType,
                            Types::Vector3<indexType> &dofs, bool set) override;

  // Paraview part
  /**
   * @brief Transfers the geometry of the element to the Paraview Adapter
   * (vtkPlotInterface).
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param paraviewAdapter[in], ParaviewAdapter object.
   * @param mainMesh [in], main mesh, normally the value 0.
   * @param subMesh[in], submesh normally the material number.
   */
  void geometryToParaview(PointerCollection &pointers,
                          vtkPlotInterface &paraviewAdapter, indexType mainMesh,
                          indexType subMesh) override;
  /**
   * @brief Compute weights for the vertices on the Paraview mesh. Required when
   * stress (or other values only available a the integration points) projection
   * is done.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param paraviewAdapter[in], ParaviewAdapter object.
   * @param mainMesh[in], main mesh, normally the value 0.
   * @param subMesh[in], submesh normally the material number.
   */
  void computeWeightsParaview(PointerCollection &pointers,
                              vtkPlotInterface &paraviewAdapter,
                              indexType mainMesh, indexType subMesh) override;
  /**
   * @brief Transfers a solution field to the Paraview Adapter
   * (vtkPlotInterface).
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param paraviewAdapter[in], ParaviewAdapter object.
   * @param mainMesh[in], main mesh, normally the value 0.
   * @param subMesh[in], submesh normally the material number.
   * @param meshId[in], mesh id of the solution field which is transferred.
   * @param order[in], order of the approximation of the solution field.
   * @param name[in], name of the solution field, predefined names are in
   * paraviewNames.
   */
  void H1SolutionToParaview(PointerCollection &pointers,
                            vtkPlotInterface &paraviewAdapter,
                            indexType mainMesh, indexType subMesh,
                            indexType meshId, indexType order,
                            std::string &name) override;
  /**
   * @brief Transfers H1 data of given order to the vertices of the Paraview
   * Adapter. Currently transfers only the linear approximation of the data.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param paraviewAdapter[in], ParaviewAdapter object.
   * @param mainMesh[in], main mesh, normally the value 0.
   * @param subMesh[in], submesh normally the material number.
   * @param Data[in], data to be transferred, Vector with the size
   * numberComponents * numberOfNodes to represent the approximation order.
   * @param numberComponents[in], number of components of the data.
   * @param order[in], order of the approximation of the data.
   * @param name[in], name of the data, predefined names are in paraviewNames.
   */
  void H1DataToParaview(PointerCollection &pointers,
                        vtkPlotInterface &paraviewAdapter, indexType mainMesh,
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
      PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
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
  void checkUpdateElement(GeometryData &geoData) override;



  /**
   * @brief Changes orientation of the face (reverts the face normal vector).
   */
  void flip() override;
  /**
   * @brief Rotates the face clockwise n times.
   *
   * @param [in] n Number of rotations.
   */
  void rotate(indexType n) override;
  /**
   * @brief Computes the mean coordinate of the face element.
   *
   * @param [in] pointers Pointer collection to global data.
   *
   * @return Mean coordinate of the face element.
   */
  auto computeMeanCoordinate(PointerCollection &pointers)
      -> Types::Vector3<prec> override;

private:
  static const GeometryTypes type;
protected:
  std::array<indexType, 4> verts, edges;
};

} // namespace HierAMuS::Geometry
