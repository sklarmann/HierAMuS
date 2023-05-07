// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "datatypes.h"
#include "pointercollection/pointercollection.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <array>
#include <forwarddeclaration.h>

#include <geometry/LinearQuadrilateral.h>
#include <types/MatrixTypes.h>

namespace HierAMuS::Geometry {

/**
 * @brief Linear Quadrilateral Geometry Element. Represents a 4-vertex face
 * element.
 *
 */
class QuadraticQuadrilateral : public LinearQuadrilateral {
public:
  QuadraticQuadrilateral();
  ~QuadraticQuadrilateral() override;

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
  auto getNumberOfVerts() -> indexType override { return 9; };
  /**
   * @brief Returns the number of edges of the element.
   * @return Number of edges as an indexType
   */
  auto getNumberOfEdges() -> indexType override { return 4; };
  /**
   * @brief Set the Vertex numbers of the object.
   *
   * @param[in] vertsIn, vector of 9 integers of type indexType refering to the
   * associated vertex numbers.
   */
  void setVerts(GeometryData &geoData, std::vector<indexType> &vertsIn) override;
  /**
   * @brief Get the global vertex numbers of the object.
   *
   * @param[out] vertsOut, vector of 9 integers of type indexType with global
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
  /**
   * @brief Prints the object to the console and to file.
   *
   * @param pointers, object containtig the pointers to global data.
   */
  void print(PointerCollection &pointers) override;

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


  /**
   * @brief Changes orientation of the face (reverts the face normal vector).
   */
  void flip() override;
  /**
   * @brief Rotates the face n times.
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
  auto getGeometryShapes(PointerCollection &pointers, IntegrationPoint &ip)
      -> H1Shapes;
  static const GeometryTypes type;
protected: 
  std::array<indexType, 5> verts_quad;
};
} // namespace HierAMuS::Geometry
