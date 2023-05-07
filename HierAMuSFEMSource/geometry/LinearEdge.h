// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <forwarddeclaration.h>
#include <geometry/Base.h>
#include <geometry/Edges.h>

#include <geometry/GeometryTypes.h>
#include <types/MatrixTypes.h>

#include <array>
#include <vector>

namespace HierAMuS::Geometry {

class LinearEdge : public Edges {
public:
  LinearEdge();
  ~LinearEdge() override;
  auto getType() -> const GeometryTypes & override;
  auto getNumberOfVerts() -> indexType override { return 2; };
  void getVerts(std::vector<indexType> &vertsOut) override;
  void getVerts(PointerCollection &pointers, std::vector<Base *> &vertsOut) override;

  auto getVerts() -> std::vector<indexType> override;

  /**
   * @brief Get the Vertex "number" of object.
   *
   * @param[in] pointers Global data collection.
   * @param[in] number Local Vertex number.
   * @return Geometry::Vertex* Pointer to vertex object.
   */
  auto getVertex(PointerCollection &pointers, indexType number) -> Vertex & override;
  /** @brief Get the global number of the vertex "number" of the object.
   *
   * @param [in] number Local vertex number.
   * @return indexType Number of Vertex.
   */
   auto getVertexNumber(indexType number) -> indexType override;
  /** @brief Check if the edge has the two vertices.
    *
    * @param[in] v1 Global number of the start vertex.
    * @param[in] v2 Global number of the end vertex.
    * @return true If the edge has the two vertices.
    * @return false If the edge does not have the two vertices.
    */
   auto hasVertices(indexType v1, indexType v2) -> bool override;
  void setVerts(GeometryData &geoData, std::vector<indexType> &vertsIn) override;
  void getNodes(PointerCollection &pointers,
                std::vector<GenericNodes *> &nodeVector, indexType meshId) override;
  /**
   * @brief Creates a std::vector of GenericNodes pointers of Internal nodes.
   *
   * @param pointers Pointer to global data.
   * @param[out] nodeVector std::vector containing the nodes of the NodeSet with
   * mesh id "meshId".
   * @param[in] meshId The mesh id of the NodeSet containing the nodes.
   */
   void getNodesInternal(PointerCollection &pointers,
                                std::vector<GenericNodes *> &nodeVector,
                                indexType meshId) override;
  void print(PointerCollection &pointers) override;

  auto getEdgeOrientation(indexType startVertex,
                          indexType endVertex) -> prec override;
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

  /**
   * @brief Computes the physical coordinates at local coordinate xi.
   *
   * @param pointers Global data collection.
   * @param xi Local coordinate xi.
   * @return Returns a vector with the x, y, z coordinates.
   */
  auto getCoordinates(PointerCollection &pointers, prec xi)
      -> Types::Vector3<prec> override;
  auto getCoordinates(PointerCollection& pointers, IntegrationPoint& IntPoint) -> Types::Vector3<prec> override;

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

  // H1Shapes
  void setH1Shapes(PointerCollection &pointers, indexType meshId,
                   indexType order, NodeTypes type) override;
  void setH1ShapesInternal(PointerCollection &pointers, indexType meshId,
                           indexType order, NodeTypes type) override;

  void getH1Dofs(PointerCollection &pointers,
                 std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                 indexType order) override;
  void getH1DofsInternal(PointerCollection &pointers,
                         std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) override;

  auto getH1Nodes(PointerCollection &pointers,
                                         indexType meshID,
                                         indexType order) -> std::vector<GenericNodes *> override;

  auto getH1NodesInternal(PointerCollection &pointers,
                                                 indexType meshID,
                                                 indexType order) -> std::vector<GenericNodes *> override;
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
  void getH1Shapes(PointerCollection &pointers, indexType order,
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
  void getH1ShapesInternal(PointerCollection &pointers, indexType order,
                           Types::VectorX<prec> &shape,
                           Types::VectorX<prec> &shapeDerivative, prec xsi) override;


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
  auto getH1Shapes(PointerCollection &pointers, indexType order,
                   IntegrationPoint &integration_point) -> H1Shapes override;
  /**
   * @brief Routine for the internal shape functions of the edge.
   *
   * @tparam prec
   * @tparam indexType
   * @param[in] pointers Collection of global data.
   * @param[in] order Order of the shape functions.
   * @param[in] integration_point Current integration point.
   * @return H1Shapes, Struct containing the shape functions and derivatives in parameter space.
   */
  auto getH1ShapesInternal(PointerCollection &pointers, indexType order,
                   IntegrationPoint &integration_point) -> H1Shapes override;

  auto getJacobian(PointerCollection &pointers, prec xi) -> prec override;
  auto getJacobian(PointerCollection &pointers,
                           IntegrationPoint &IntegrationPt)
      -> Types::MatrixXX<prec> override;

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
                              ShapeFunctionTypes shapeType,
                              indexType shapeOrder,
                              Types::VectorX<prec> &Solution,
                              indexType propNumber,
                              Types::VectorX<prec> &direction, bool local,
                              bool add) override;


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
                           vtkPlotInterface &paraviewAdapter,
                           indexType mainMesh, indexType subMesh) override;
   /**
    * @brief Compute weights for the vertices on the Paraview mesh. Required
    * when stress (or other values only available a the integration points)
    * projection is done.
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


  auto getDirectionVector(PointerCollection &pointers) -> Types::Vector3<prec> override;

  auto getA1Vector(PointerCollection &pointers,
                   IntegrationPoint &integration_point) -> Types::Vector3<prec> override;

  void setAllNodeBoundaryConditionMeshId(PointerCollection &pointers,
                                                 indexType meshId,
                                                 indexType dof) override;

  /**
   * @brief Changes orientation of the face (reverts the face normal vector).
   */
  void flip() override;

private:
  std::array<indexType, 2> m_verts;
  static const GeometryTypes type;
};

} // namespace HierAMuS
