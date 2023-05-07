// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "geometry/Edges.h"
#include "geometry/Vertex.h"
#include "pointercollection/pointercollection.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <forwarddeclaration.h>
#include <geometry/Base.h>

#include <geometry/GeometryTypes.h>
#include <tuple>
#include <types/MatrixTypes.h>

#include <vector>

namespace HierAMuS::Geometry {
/**
 * @brief Face Orientations
 *
 * This enum represents the face orientation flags.
 */
enum class faceorientation {
  p_1, // positive normal direction and edge vertices match with first face edge
  p_2, // positive normal direction and rotated by once
  p_3, // positive normal direction and rotated by twice
  p_4, // positive normal direction and rotated by three times
  n_1, // negative normal direction and edge vertices match with first face edge
  n_2, // negative normal direction and rotated once
  n_3, // negative normal direction and rotated twice
  n_4  // negative normal direction and rotated three times
};
class Faces : public Base {
public:
  Faces();
  ~Faces() override;
  auto getGroupType() -> const GeometryTypes & override;

  virtual auto getCoordinates(PointerCollection &pointers,
                              IntegrationPoint &integrationPoint)
      -> Types::Vector3<prec> override {
    return {};
  };


  /**
   * @brief Returns the number of edges of the element. 
   * @return Number of edges as an indexType
   */
  virtual auto getNumberOfEdges() -> indexType = 0;

  virtual auto getTangent_G1(PointerCollection &pointers,
                             IntegrationPoint &integrationPoint)
      -> Types::Vector3<prec> {
    return {};
  };

  virtual auto getTangent_G2(PointerCollection &pointers,
                             IntegrationPoint &integrationPoint)
      -> Types::Vector3<prec> {
    return {};
  };

  virtual auto getVertex(PointerCollection &pointers, indexType local_number)
      -> Geometry::Vertex * {
    return nullptr;
  };
  virtual auto getEdge(PointerCollection &pointers, indexType local_number)
      -> Geometry::Edges * {
    return nullptr;
  };

  /** @brief Returns the global edge number based on the local edge number of the face.
   *
   *  detailed description
   *
   *  @param[in] local_number local edge number of the face
   *  @return indexType global edge number
   */
  virtual auto getEdge(indexType local_number) -> indexType = 0;

  /** @brief Check if the face has the three vertices.
    *
    * @param[in] v1 Global number of a corner vertex.
    * @param[in] v2 Global number of a corner vertex.
    * @param[in] v3 Global number of a corner vertex.
    * @return true If the face has the three vertices.
    * @return false If the face does not have the three vertices.
    */
  virtual auto hasVertices(indexType v1, indexType v2, indexType v3) -> bool = 0;

  virtual auto getOrientation(PointerCollection &pointers, indexType vertex1,
                              indexType vertex2) -> faceorientation {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling getOrientation for Faces, not implemented!");
    return {};
  };

  virtual void modifyIntegrationpoint(IntegrationPoint &IP, prec &shapeFactor,
                                      faceorientation orientation) {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling modifyInegrationpoint for Faces, not implemented!");
  };

  // H1 Shapes
  virtual void setH1Shapes(PointerCollection &pointers, indexType meshId,
                           indexType order, NodeTypes type) {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling setH1Shapes for Faces, not implemented!");
  };
  virtual void setH1ShapesInternal(PointerCollection &pointers,
                                   indexType meshId, indexType order,
                                   NodeTypes type) {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling setH1ShapesInternal for Faces, not implemented!");
  };

  virtual auto getH1Dofs(PointerCollection &pointers, indexType meshID,
                         indexType order) -> std::vector<DegreeOfFreedom *> {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling getH1Dofs for Faces, not implemented!");
    return {};
  };

  virtual void getH1Dofs(PointerCollection &pointers,
                         std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling getH1Dofs for Faces, not implemented!");
  };

  virtual void getH1DofsInternal(PointerCollection &pointers,
                                 std::vector<DegreeOfFreedom *> &Dofs,
                                 indexType meshID, indexType order) {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling getH1DofsInternal for Faces, not implemented!");
  };

  virtual auto getH1Nodes(PointerCollection &pointers, indexType meshID,
                          indexType order) -> std::vector<GenericNodes *> {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling getH1Nodes for Edges, not implemented!");
    return {};
  };

  virtual auto getH1NodesInternal(PointerCollection &pointers, indexType meshID,
                                  indexType order)
      -> std::vector<GenericNodes *> {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling getH1Nodes for Edges, not implemented!");
    return {};
  };

  virtual auto getHDivNodes(PointerCollection &pointers, indexType meshID,
                            indexType order) -> std::vector<GenericNodes *> {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling getH1Nodes for Edges, not implemented!");
    return {};
  };

  virtual auto getHDivNodesInternal(PointerCollection &pointers,
                                    indexType meshID, indexType order)
      -> std::vector<GenericNodes *> {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling getH1Nodes for Edges, not implemented!");
    return {};
  };

  virtual auto getH1Shapes(PointerCollection &pointers, indexType order,
                           IntegrationPoint &IntegrationPt) -> H1Shapes {
    HierAMuS::Geometry::Faces::throwError(
        "Function getH1Shapes for the geometry element not implemented!");
    return {};
  };

  virtual auto
  getH1ShapesInternal(PointerCollection &pointers, indexType order,
                      IntegrationPoint &IntegrationPt,
                      faceorientation orientation = faceorientation::p_1)
      -> H1Shapes {
    HierAMuS::Geometry::Faces::throwError(
        "Function getH1Shapes for the geometry element not implemented!");
    return {};
  };

  virtual void getH1Shapes(PointerCollection &pointers, indexType order,
                           Types::VectorX<prec> &shape,
                           Types::Matrix2X<prec> &shapeDerivative, prec xsi,
                           prec eta) {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling getH1Shapes for Faces, not implemented!");
  };

  virtual void getH1ShapesInternal(PointerCollection &pointers, indexType order,
                                   Types::VectorX<prec> &shape,
                                   Types::Matrix2X<prec> &shapeDerivative,
                                   prec xsi, prec eta) {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling getH1ShapesInternal for Faces, not implemented!");
  };

  // HDivShapes
  virtual void setHDivShapes(PointerCollection &pointers, indexType meshid,
                             indexType order, NodeTypes type) {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling setHDivShapes for Faces, not implemented!");
  };

  virtual void getHDivDofs(PointerCollection &pointers,
                           std::vector<DegreeOfFreedom *> &Dofs,
                           indexType meshID, indexType order, NodeTypes type) {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling getHDivDofs for Faces, not implemented!");
  };
  virtual void getHDivShapes(PointerCollection &pointers,
                             indexType order,
                             Types::Matrix2X<prec> &shape,
                             Types::VectorX<prec> &dshape, prec xi,
                             prec eta) {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling getHDivShapes for Faces, not implemented!");
  };

  virtual auto getHDivShapes(PointerCollection &pointers, indexType order,
                             IntegrationPoint &IntegrationPt) -> HDivShapes {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling getHDivShapes for Faces, not implemented!");
    return {};
  };

  // L2 shapes

  virtual void getL2Dofs(PointerCollection &pointers,
                         std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order){};

  virtual auto getL2Shapes(PointerCollection &pointers, indexType order,
                           IntegrationPoint &IntegrationPt) -> L2Shapes {
    return {};
  };

  virtual void getL2ShapesInternal(PointerCollection &pointers, indexType order,
                                   Types::VectorX<prec> &shape,
                                   Types::Matrix2X<prec> &shapeDerivative,
                                   prec xsi, prec eta){};

  virtual void setL2Shapes(PointerCollection &pointers, indexType meshId,
                           indexType order, NodeTypes type){};

  // Special Plate Shapes
  virtual void setSpecialPlateShapes(PointerCollection &pointers,
                                     indexType meshid, indexType order,
                                     NodeTypes type) {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling setSpecialPlateShapes for Faces, not implemented!");
  };
  virtual auto getSpecialPlateDofs(PointerCollection &pointers,
                                   indexType meshID, indexType order,
                                   NodeTypes type)
      -> std::vector<DegreeOfFreedom *> {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling setSpecialPlateShapes for Faces, not implemented!");
    return {};
  };
  virtual auto getSpecialPlateShapes(PointerCollection &pointers,
                                     IntegrationPoint &intPoint,
                                     indexType order) -> SpecialPlateShapes {
    HierAMuS::Geometry::Faces::throwError(
        "Error when calling getSpecialPlateShapes for Faces, not implemented!");
    return {};
  };

  auto getJacobian(PointerCollection &pointers, IntegrationPoint &IntegrationPt)
      -> Types::MatrixXX<prec> override {
    throw std::runtime_error(
        "Function getJacobian not implemented for the used geometry object!");
    return {};
  }

  void setAllNodeBoundaryConditionMeshId(PointerCollection &pointers,
                                         indexType meshId,
                                         indexType dof) override;

  // Paraview
  void geometryToParaview(PointerCollection &pointers,
                          vtkPlotInterface &paraviewAdapter, indexType mainMesh,
                          indexType subMesh) override{};
  void computeWeightsParaview(PointerCollection &pointers,
                              vtkPlotInterface &paraviewAdapter,
                              indexType mainMesh, indexType subMesh) override{};
  void H1SolutionToParaview(PointerCollection &pointers,
                            vtkPlotInterface &paraviewAdapter,
                            indexType mainMesh, indexType subMesh,
                            indexType meshId, indexType order,
                            std::string &name) override{};
  void H1DataToParaview(PointerCollection &pointers,
                        vtkPlotInterface &paraviewAdapter, indexType mainMesh,
                        indexType subMesh, Types::VectorX<prec> &Data,
                        indexType numberComponents, indexType order,
                        std::string &name) override{};

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
  virtual void checkUpdateElement(GeometryData &geoData) = 0;

  /**
   * @brief Changes orientation of the face (reverts the face normal vector).
   */
  virtual void flip() = 0;
  /**
   * @brief Rotates the face n times.
   *
   * @param [in] n Number of rotations.
   */
  virtual void rotate(indexType n) = 0;
  /**
   * @brief Computes the mean coordinate of the face element.
   *
   * @param [in] pointers Pointer collection to global data.
   *
   * @return Mean coordinate of the face element.
   */
  virtual auto computeMeanCoordinate(PointerCollection &pointers)
      -> Types::Vector3<prec> = 0;

private:
  static const GeometryTypes type;
  static void throwError(const std::string &msg) {
    throw std::runtime_error(msg);
  };
};

} // namespace HierAMuS::Geometry
