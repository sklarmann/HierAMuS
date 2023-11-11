// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "geometry/VertexRuntime.h"
#include "geometry/Volumes/VolumesH1Interface.h"
#include "geometry/Volumes/VolumesRuntimeDataInterface.h"
#include <geometry/GeometryTypes.h>

#include <types/MatrixTypes.h>

#include <vector>

namespace HierAMuS::Geometry {
class LinearBrickData;
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
class LinearBrickRuntime
    : public VolumesRuntimeDataInterface<8, 12, 6, LinearBrickData,
                                         LinearBrickRuntime>,
      public VolumesH1Interface {

public:
  LinearBrickRuntime(GeometryData &geoData, LinearBrickData &data_element);
  ~LinearBrickRuntime() override;

  auto getH1Volume() -> VolumesH1Interface * override { return this; };

  auto getType() -> const GeometryTypes & override;
  void print(spdlog::logger &Log) override;

  /**
   * @brief Get the global edge numbers of the geometric element.
   *
   * @param edgesOut[out], a vector with the global numbers of the Edges of type
   * indexType and of size 12.
   */
  void getEdgeNumbers(std::vector<indexType> &edgesOut) override;

  /**
   * @brief Get the global face numbers of the geometric element.
   *
   * @param facesOut[out], a vector with the global numbers of the Faces of type
   * indexType and of size 6.
   */
  void getFaceNumbers(std::vector<indexType> &facesOut) override;

  /**
   * @brief Get the Integration Points for the brick element.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @return IntegrationPoints
   */
  auto getIntegrationPoints(indexType elementId) -> IntegrationPoints override;

  // Geometric mapping

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

  void getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                 indexType order) override;
  void getH1DofsInternal(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) override;

  /**
   * @brief Returns a vector of the H1 nodes for the given meshId and shape
   * function order.
   *
   * @param pointers, object containing the pointers to global data.
   * @param meshID, mesh id of the nodes for the H1 shape functions.
   * @param order, order of the H1 shape functions.
   * @return std::vector<GenericNodes *>, vector of the GenericNodes.
   */
  auto getH1Nodes(indexType meshID, indexType order)
      -> std::vector<GenericNodes *> override;
  auto getH1NodesInternal(indexType meshID, indexType order)
      -> std::vector<GenericNodes *> override;

  void getH1Shapes(indexType order,
                   Types::VectorX<prec> &shape,
                   Types::Matrix3X<prec> &shapeDerivative, prec xsi, prec eta,
                   prec zeta) override;

  void getH1ShapesInternal(indexType order,
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
  auto getH1Shapes(indexType order,
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
  auto getH1ShapesInternal(indexType order,
                           IntegrationPoint &IntegrationPt)
      -> H1Shapes override;

  // Paraview
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
   * @param paraviewAdapter[in], ParaviewAdapter object.
   * @param mainMesh[in], main mesh, normally the value 0.
   * @param subMesh[in], sub mesh normally the material number.
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

private:
  LinearBrickData &m_LinearBrick_Data;

  static const GeometryTypes type;
};

} // namespace HierAMuS::Geometry