// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once



#include <forwarddeclaration.h>

#include <equations/GenericNodes.h>
#include <equations/Nodetypes.h>
#include <geometry/GeometryTypes.h>
#include <shapefunctions/IntegrationsPoints/IntegrationPoints.h>
#include <stdexcept>
#include <types/MatrixTypes.h>
#include <vector>

#include <tuple>

namespace HierAMuS::Geometry {

struct H1Shapes {
  Types::VectorX<prec> shapes;
  Types::MatrixXX<prec> shapeDeriv;
};

struct L2Shapes {
    Types::VectorX<prec> shapes;
};


struct HDivShapes {
  Types::MatrixXX<prec> shapes;
  Types::MatrixXX<prec> shapeDeriv;
};

struct SpecialPlateShapes {
  Types::Matrix2X<prec> shapes;
  Types::Matrix2X<prec> shapeDx;
  Types::Matrix2X<prec> shapeDy;
};


class GeometryData;

/**
 * @brief Base geometry class.
 * All other geometry classes are derived from this class.
 * This class handles the connection to the degrees of freedom of a geometric
 * object.
 *
 */
class Base {
  using ptrCol = PointerCollection;

public:
  Base();
  virtual ~Base();

  /**
   * @brief Set the id of the geometric object for identification.
   *
   * @param[in] id Id of the object.
   */
  void setId(indexType id) { this->id = id; };
  /**
   * @brief Get the id of the current object.
   *
   * @return Returns the if of the object of type indextype.
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
   * @param[in] set If true, boudary contidions are overridden, otherwise only 1
   * is considered.
   */
  virtual void setBoundaryCondition(PointerCollection &pointers,
                                    indexType meshId, indexType order,
                                    ShapeFunctionTypes shapeType,
                                    Types::Vector3<indexType> &dofs,
                                    bool set){};

  /**
   * @brief Set a NodeSet for the geometric element containing the nodes.
   *
   * @param pointers Global data collection.
   * @param[in] meshID The mesh id of the nodes.
   * @param[in] numberOfNodes Number of nodes the NodeSet contains.
   * @param[in] type The type of the nodes.
   */
  void setNodeSet(PointerCollection &pointers, indexType meshID,
                  indexType numberOfNodes, NodeTypes type);

  /**
   * @brief Returns a pointer to the NodeSet "setNumber" with zero based
   * indexing.
   *
   * @param pointers Pointers to global data.
   * @param[in] setNumber The set number with zero based indexing.
   * @return The requested NodeSet.
   */
  auto getSet(ptrCol &pointers, indexType setNumber) -> NodeSet *;
  /**
   * @brief Returns a pointer to the node set with mesh id "meshId".
   *
   * @param pointers Pointers to global data.
   * @param[in] meshId The mesh id of the set requested.
   * @return The requested node set with mesh id "meshId".
   */
  auto getSetMeshId(ptrCol &pointers, indexType meshId) -> NodeSet *;
  /**
   * @brief Creates a vector with pointers to all NodeSet the geometric element
   * has.
   *
   * @param pointers Pointers to global data.
   * @param[out] sets std::vector of NodeSet pointers.
   */
  void getSets(ptrCol &pointers, std::vector<NodeSet *> &sets);
  /**
   * @brief Creates a std::vector of GenericNodes pointers.
   *
   * @param pointers Pointer to global data.
   * @param[out] nodeVector std::vector containing the nodes of the NodeSet with
   * mesh id "meshId".
   * @param[in] meshId The mesh id of the NodeSet containing the nodes.
   */
  void getNodesOfSet(ptrCol &pointers, std::vector<GenericNodes *> &nodeVector,
                     indexType meshId);
  auto getNodesOfSet(ptrCol &pointers, indexType meshId) -> std::vector<GenericNodes*>;
  /**
   * @brief Creates a std::vector of GenericNodes pointers.
   *
   * @param pointers Pointer to global data.
   * @param[out] nodeVector std::vector containing the nodes of the NodeSet with
   * mesh id "meshId".
   * @param[in] meshId The mesh id of the NodeSet containing the nodes.
   */
  virtual void getNodes(ptrCol &pointers,
                        std::vector<GenericNodes *> &nodeVector,
                        indexType meshId);
  /**
   * @brief Creates a std::vector of GenericNodes pointers of Internal nodes.
   *
   * @param pointers Pointer to global data.
   * @param[out] nodeVector std::vector containing the nodes of the NodeSet with
   * mesh id "meshId".
   * @param[in] meshId The mesh id of the NodeSet containing the nodes.
   */
  virtual void getNodesInternal(ptrCol &pointers,
                                std::vector<GenericNodes *> &nodeVector,
                                indexType meshId);
  /**
   * @brief Check if a set with the given meshId is present.
   *
   * @param[in] pointers Pointers to global data.
   * @param[in] meshId Mesh id to check.
   * @return Returns true if a set is present, otherwise false.
   */
  auto checkSetWithMeshId(ptrCol &pointers, indexType meshId) -> bool;

  /**
   * @brief Sets boundary condition at degree of freedom "dof" on all nodes in
   * NodeSet with mesh id "meshid".
   *
   * @param pointers Pointers to global data.
   * @param meshId Mesh id of the NodeSet.
   * @param dof Degree of freedom number to set the boundary condition, zero
   * based indexing.
   */
  virtual void setAllNodeBoundaryConditionMeshId(ptrCol &pointers,
                                                 indexType meshId,
                                                 indexType dof);
  void getAllEquationsIds(ptrCol &pointers,
                          std::vector<DegreeOfFreedom *> &Dofs);
  void getNodeEquationIds(ptrCol &pointers,
                          std::vector<DegreeOfFreedom *> &Dofs,
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
  virtual void setLoad(ptrCol &pointers, indexType meshid,
                       ShapeFunctionTypes shapeType, indexType shapeOrder,
                       Types::VectorX<prec> &Loads, indexType propNumber,
                       Types::VectorX<prec> &direction, bool local, bool add){};

  /**
   * @brief Sets the solution using a Geometric Object.
   * If boundary conditions are not set, they will be set, if the prescribed value is not zero.
   *
   * @param[in] pointers Pointers to global data.
   * @param[in] meshid Mesh id of the degrees of freedom to add the load.
   * @param[in] shapeType Shape function Type (H0 or H1).
   * @param[in] shapeOrder Order of the shape functions to approximate delta u.
   * @param[in] Solution Vector with prescribed solution, length must be a multiple of 3,
   * depending on the length, the load interpolation function is chosen.
   * @param[in] propNumber Proportional load function number for the load.
   * @param[in] direction If local is a function, this vector estimates if the
   * integration is in positive or negative direction of an edge.
   * @param[in] local If true, the loads are treated in the local coordinate
   * system.
   * @param[in] add If true, the loads will be added to the current prescribed values
   * otherwise prescribed values will be overridden.
   */
  virtual void setPrescribedSolution(ptrCol &pointers, indexType meshid,
                       ShapeFunctionTypes shapeType, indexType shapeOrder,
                       Types::VectorX<prec> &Solution, indexType propNumber,
                       Types::VectorX<prec> &direction, bool local, bool add){};

  /**
   * @brief Returns the type of the geometric object.
   *
   * @return Type of the object.
   */
  virtual auto getType() -> const GeometryTypes &;
  /**
   * @brief Returns the group typ of the geometric element.
   * The group types are: Edges, Faces, Volumes.
   *
   * @return GeometryTypes as the type
   */
  virtual auto getGroupType() -> const GeometryTypes &;

  // Output Information
  virtual void print(PointerCollection &pointers) = 0;
  void printEqInfo(PointerCollection &pointers);

  virtual void setCoordinates(prec x = 0, prec y = 0, prec z = 0){};
  virtual void setCoordinates(const Types::Vector3<prec> &coorIn){};

  virtual auto getCoordinates() -> Types::Vector3<prec>
  {
    return {};
  };
  virtual auto getCoordinates(PointerCollection& pointers,
                              prec xi) -> Types::Vector3<prec>
  {
    return {};
  };
  virtual auto getCoordinates(PointerCollection& pointers,
                              prec xi, prec eta) -> Types::Vector3<prec>
  {
    return {};
  };
  virtual auto getCoordinates(PointerCollection& pointers,
                              prec xi, prec eta, prec zeta) -> Types::Vector3<prec>
  {
    return {};
  };


  virtual auto getCoordinates(PointerCollection &pointers,
                              IntegrationPoint &IntPoint) -> Types::Vector3<prec> = 0;

  /**
   * @brief Get the Face Normal vector
   * Only virtual method, must be implemented in the derived face classes.
   *
   * @param pointers[in], pointers to global data.
   * @return Types::Vector3<prec>, Normal vector of the face
   */
  virtual auto getFaceNormal(PointerCollection &pointers) -> Types::Vector3<prec>{
    throw std::runtime_error(
        "Function getFaceNormal not implemented for the used geometry object!");
        return {};
  };

  // Getters and setters
  virtual void setVerts(GeometryData &geoData, std::vector<indexType> &vertsIn);
  virtual auto getNumberOfVerts() -> indexType { return 0; };
  virtual void getVerts(std::vector<indexType> &vertsOut) {
    vertsOut.resize(0);
  };
  virtual void getVerts(ptrCol &pointers, std::vector<Base *> &vertsOut) {
    vertsOut.resize(0);
  };

  virtual auto getVerts() -> std::vector<indexType> { return {}; };

  virtual void setEdges(const std::vector<indexType> &edgesIn){};
  virtual void getEdges(std::vector<indexType> &edgesOut) {
    edgesOut.resize(0);
  };
  virtual void getEdges(ptrCol &pointers, std::vector<Base *> &edgesOut) {
    edgesOut.resize(0);
  };

  virtual void setFaces(std::vector<indexType> &facesIn){};
  virtual auto getNumberOfFaces() -> indexType { return 0; };
  virtual void getFaces(std::vector<indexType> &facesOut) {
    facesOut.resize(0);
  };
  virtual void getFaces(ptrCol &pointers, std::vector<Base *> &facesOut) {
    facesOut.resize(0);
  };

  virtual auto getIntegrationPoints(ptrCol &pointers, indexType elementId) -> IntegrationPoints;

  virtual auto getJacobian(ptrCol& pointers,
                           IntegrationPoint& IntegrationPt) -> Types::MatrixXX<prec>
  {
    throw std::runtime_error(
        "Function getJacobian not implemented for the used geometry object!");
    return {};
  }

  virtual void geometryToParaview(ptrCol &pointers,
                                  vtkPlotInterface &paraviewAdapter,
                                  indexType mainMesh, indexType subMesh){};
  virtual void computeWeightsParaview(ptrCol &pointers,
                                      vtkPlotInterface &paraviewAdapter,
                                      indexType mainMesh, indexType subMesh){};
  virtual void H1SolutionToParaview(ptrCol &pointers,
                                    vtkPlotInterface &paraviewAdapter,
                                    indexType mainMesh, indexType subMesh,
                                    indexType meshId, indexType order,
                                    std::string &name){};
  virtual void H1DataToParaview(PointerCollection &pointers,
                                vtkPlotInterface &paraviewAdapter,
                                indexType mainMesh, indexType subMesh,
                                Types::VectorX<prec> &Data,
                                indexType numberComponents, indexType order,
                                std::string &name){};
  virtual void projectDataToParaviewVertices(PointerCollection &pointers,
                                          vtkPlotInterface &paraviewAdapter,
                                          indexType mainMesh, indexType subMesh,
                                          indexType order,
                                          IntegrationPoint &IntegrationPt,
                                          Types::VectorX<prec> &data,
                                          indexType numberComponents,
                                          std::string name){};

protected:
  indexType id;
  unsigned char numberOfNodeSets;
  indexType nodeSetId;
  bool nodeSetInit;

private:
  static const HierAMuS::Geometry::GeometryTypes type;
};

} // namespace HierAMuS
