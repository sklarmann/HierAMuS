// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "pointercollection/pointercollection.h"
#include <forwarddeclaration.h>
#include <geometry/GeometryTypes.h>
#include <geometry/Volumes.h>

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
class LinearBrick : public Volumes {
public:
  LinearBrick();
  ~LinearBrick() override;
  auto getType() -> const GeometryTypes & override;
  void print(PointerCollection &pointers) override;

  /**
   * @brief Set the Vertex numbers of the geometric element.
   * They need to be given according to the numbering in the picture above.
   *
   * @param vertsIn[in], a vector with the global numbers of the vertices of
   * type indexType and of size 8.
   */
  void setVerts(GeometryData &geoData, std::vector<indexType> &vertsIn) override;
  /**
   * @brief Set the Edge numbers of the geometric element.
   * They need to be given according to the numbering in the picture above.
   *
   * @param vertsIn[in], a vector with the global numbers of the Edges of type
   * indexType and of size 12.
   */
  void setEdges(const std::vector<indexType> &edgesIn) override;
  /**
   * @brief Set the Face numbers of the geometric element.
   * They need to be given according to the numbering in the picture above.
   *
   * @param vertsIn[in], a vector with the global numbers of the Faces of type
   * indexType and of size 6.
   */
  void setFaces(std::vector<indexType> &facesIn) override;

  /**
   * @brief Get the global vertex numbers of the geometric element.
   *
   * @param vertsOut[out], a vector with the global numbers of the vertices of
   * type indexType and of size 8.
   */
  void getVerts(std::vector<indexType> &vertsOut) override;

  void getVerts(PointerCollection &pointers,
                std::vector<Base *> &vertsOut) override;

  /**
   * @brief Get coordinate at integration point.
   *
   * @param pointers, pointer to global data.
   * @param IntPoint, integration point to evaluate coordiantes at.
   *
   * @return Vector3 with the x,y,z coordinates.
   */
  auto getCoordinates(PointerCollection& pointers, IntegrationPoint& IntPoint) -> Types::Vector3<prec> override;
  /**
   * @brief Get the Number Of Verts of the geometric object, return 8
   *
   * @return indexType, returns 8.
   */
  auto getNumberOfVerts() -> indexType override { return 8; };

  /**
   * @brief Get the global edge numbers of the geometric element.
   *
   * @param edgesOut[out], a vector with the global numbers of the Edges of type
   * indexType and of size 12.
   */
  void getEdges(std::vector<indexType> &edgesOut) override;
  /**
   * @brief Get a vector with pointers of type Geometry::Base to the edges of
   * the geometric element.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param edgesOut[out], a vector of size 12 with pointers to the edges of the
   * geometric element.
   */
  void getEdges(PointerCollection &pointers,
                std::vector<Base *> &edgesOut) override;
  /**
   * @brief Get the Number of Edges of the element, returns 12.
   *
   * @return indexType, returns 12.
   */
  auto getNumberOfEdges() -> indexType override { return 12; };
  /**
   * @brief Get the global face numbers of the geometric element.
   *
   * @param facesOut[out], a vector with the global numbers of the Faces of type
   * indexType and of size 6.
   */
  void getFaces(std::vector<indexType> &facesOut) override;
  /**
   * @brief Get a vector with pointers of type Geometry::Base to the faces of
   * the geometric element.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param facesOut[out], a vector of size 6 with pointers to the faces of the
   * geometric element.
   */
  void getFaces(PointerCollection &pointers,
                std::vector<Base *> &facesOut) override;

  /**
   * @brief Get the Integration Points for the brick element.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @return IntegrationPoints
   */
  auto getIntegrationPoints(PointerCollection &pointers, indexType elementId)
  -> IntegrationPoints override;

  // Geometric mapping
  void getJacobian(PointerCollection &pointers, Types::Matrix33<prec> &jacobi,
                   prec xsi, prec eta, prec zeta) override;
  /**
   * @brief Get the Jacobian matrix at the current IntegrationPoint. New version
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param IntegrationPt[in], current IntegrationPoint.
   * @return Types::MatrixXX<prec>, Jacobian matrix at the current
   * IntegrationPoint.
   */
  auto getJacobian(PointerCollection &pointers, IntegrationPoint &IntegrationPt)
      -> Types::MatrixXX<prec> override;

  // H1 Shapes
  /**
   * @brief Assignes the nodes to the geometry object for H1 shape functions
   * with given order associated with the specific meshId. A single element can
   * have multiple H1 nodes with different meshids (e.g. for different solutions
   * like displacements and rotations).
   *
   * @param pointers, object containing the pointers to global data.
   * @param meshId, mesh id of the nodes for the H1 shape functions.
   * @param order, order of the H1 shape functions.
   * @param type, type of the degrees of freedom, currently unused.
   */
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

  /**
   * @brief Returns a vector of the H1 nodes for the given meshId and shape
   * function order.
   *
   * @param pointers, object containtig the pointers to global data.
   * @param meshID, mesh id of the nodes for the H1 shape functions.
   * @param order, order of the H1 shape functions.
   * @return std::vector<GenericNodes *>, vector of the GenericNodes.
   */
  auto getH1Nodes(PointerCollection &pointers, indexType meshID,
                  indexType order) -> std::vector<GenericNodes *> override;
  auto getH1NodesInternal(PointerCollection &pointers, indexType meshID,
                          indexType order)
      -> std::vector<GenericNodes *> override;

  void getH1Shapes(PointerCollection &pointers, indexType order,
                   Types::VectorX<prec> &shape,
                   Types::Matrix3X<prec> &shapeDerivative, prec xsi, prec eta,
                   prec zeta) override;

  void getH1ShapesInternal(PointerCollection &pointers, indexType order,
                           Types::VectorX<prec> &shape,
                           Types::Matrix3X<prec> &shapeDerivative, prec xsi,
                           prec eta, prec zeta) override;

  /**
   * @brief New version of getH1Shapes.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param order[in], order of the H1 shape functions.
   * @param IntegrationPt[in], integration point.
   * @return H1Shapes, H1 shape functions.
   */
  auto getH1Shapes(PointerCollection &pointers, indexType order,
                   IntegrationPoint &IntegrationPt) -> H1Shapes override;
  /**
   * @brief New version of getH1ShapesInternal.
   * Computes and returns the H1 bubble functions of the volume element.
   *
   * @todo Implementation of higher order H1 shapes.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param order[in], order of the H1 shape functions.
   * @param IntegrationPt[in], integration point.
   * @return H1Shapes, H1 shape functions.
   */
  auto getH1ShapesInternal(PointerCollection &pointers, indexType order,
                           IntegrationPoint &IntegrationPt)
      -> H1Shapes override;

  // Paraview
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
   * @brief Projects the the data at the current Integration point to the vertices in the Paraview mesh.
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


  void setAllNodeBoundaryConditionMeshId(PointerCollection &pointers,
                                     indexType meshId, indexType dof) override;

private:
  std::array<indexType, 8> m_verts;
  std::array<indexType, 12> m_edges;
  std::array<indexType, 6> m_faces;
  static const GeometryTypes type;
};

} // namespace HierAMuS