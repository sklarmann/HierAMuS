// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MeshIdNodeList.h"
#include "geometry/Edges/EdgesData.h"
#include "geometry/VertexData.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

#include <geometry/GeometryBaseData.h>

#include <geometry/GeometryTypes.h>
#include <tuple>
#include <types/MatrixTypes.h>

#include <vector>
#include "FaceOrientationFlags.h"

namespace HierAMuS {
class EquationHandler;
}

namespace HierAMuS::Geometry {
class FacesRuntime;

class FacesData : public GeometryBaseData {
public:
  FacesData();
  ~FacesData() override;
  auto getGroupType() -> const GeometryTypes & override;

  virtual auto getCoordinates(IntegrationPoint &integrationPoint)
      -> Types::Vector3<prec> = 0;

  virtual auto getCoordinates(prec xi, prec eta)
      -> Types::Vector3<prec> = 0;

  virtual auto getRuntimeObject(GeometryData &geoData) -> std::shared_ptr<FacesRuntime> = 0;
  /**
   * @brief Returns the number of edges of the element.
   * @return Number of edges as an indexType
   */
  virtual auto getNumberOfEdges() -> indexType = 0;

  virtual auto getTangent_G1(IntegrationPoint &integrationPoint)
      -> Types::Vector3<prec> {
    return {};
  };

  virtual auto getTangent_G2(IntegrationPoint &integrationPoint)
      -> Types::Vector3<prec> {
    return {};
  };

  virtual auto getVertex(indexType local_number)
      -> Geometry::VertexData * = 0;


  /** @brief Returns the global edge number based on the local edge number of
   * the face.
   *
   *  detailed description
   *
   *  @param[in] local_number local edge number of the face
   *  @return indexType global edge number
   */
  virtual auto getEdge(indexType local_number) -> EdgesData * = 0;

  /** @brief Check if the face has the three vertices.
   *
   * @param[in] v1 Global number of a corner vertex.
   * @param[in] v2 Global number of a corner vertex.
   * @param[in] v3 Global number of a corner vertex.
   * @return true If the face has the three vertices.
   * @return false If the face does not have the three vertices.
   */
  virtual auto hasVertices(indexType v1, indexType v2, indexType v3)
      -> bool = 0;

  virtual auto getOrientation(indexType vertex1,
                              indexType vertex2) -> faceorientation {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling getOrientation for Faces, not implemented!");
    return {};
  };

  virtual void modifyIntegrationpoint(IntegrationPoint &IP, prec &shapeFactor,
                                      faceorientation orientation) {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling modifyInegrationpoint for Faces, not implemented!");
  };

  // H1 Shapes
  virtual void setH1Shapes(indexType meshId,
                           indexType order, NodeTypes type) = 0;
  virtual void setH1ShapesInternal(indexType meshId, indexType order,
                                   NodeTypes type) = 0;

  virtual auto getH1Dofs(indexType meshID,
                         indexType order) -> std::vector<DegreeOfFreedom *> {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling getH1Dofs for Faces, not implemented!");
    return {};
  };

  virtual void getH1Dofs(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling getH1Dofs for Faces, not implemented!");
  };

  virtual void getH1DofsInternal(std::vector<DegreeOfFreedom *> &Dofs,
                                 indexType meshID, indexType order) = 0;

  virtual auto getH1NodesList(indexType meshID,
                              indexType order) -> MeshIdNodeList = 0;

  virtual auto getH1Nodes(indexType meshID,
                          indexType order) -> std::vector<GenericNodes *> {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling getH1Nodes for Face, not implemented!");
    return {};
  };

  virtual auto getH1NodesInternal(indexType meshID,
                                  indexType order)
      -> std::vector<GenericNodes *> {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling getH1NodesInternal for Face, not implemented!");
    return {};
  };

  virtual auto getHDivNodes(indexType meshID,
                            indexType order) -> std::vector<GenericNodes *> {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling getH1Nodes for Edges, not implemented!");
    return {};
  };

  virtual auto getHDivNodesInternal(indexType meshID, indexType order)
      -> std::vector<GenericNodes *> {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling getH1Nodes for Edges, not implemented!");
    return {};
  };

  virtual auto getH1Shapes(indexType order,
                           IntegrationPoint &IntegrationPt) -> H1Shapes {
    HierAMuS::Geometry::FacesData::throwError(
        "Function getH1Shapes for the geometry element not implemented!");
    return H1Shapes(0, 0);
  };

  virtual auto
  getH1ShapesInternal(indexType order,
                      IntegrationPoint &IntegrationPt,
                      faceorientation orientation = faceorientation::p_1)
      -> H1Shapes {
    HierAMuS::Geometry::FacesData::throwError(
        "Function getH1Shapes for the geometry element not implemented!");
    return H1Shapes(0, 0);
  };

  virtual void getH1Shapes(indexType order,
                           Types::VectorX<prec> &shape,
                           Types::Matrix2X<prec> &shapeDerivative, prec xsi,
                           prec eta) {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling getH1Shapes for Faces, not implemented!");
  };

  virtual void getH1ShapesInternal(indexType order,
                                   Types::VectorX<prec> &shape,
                                   Types::Matrix2X<prec> &shapeDerivative,
                                   prec xsi, prec eta) {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling getH1ShapesInternal for Faces, not implemented!");
  };

  // HDivShapes
  virtual void setHDivShapes(indexType meshid,
                             indexType order, NodeTypes type) {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling setHDivShapes for Faces, not implemented!");
  };

  virtual void getHDivDofs(std::vector<DegreeOfFreedom *> &Dofs,
                           indexType meshID, indexType order, NodeTypes type) {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling getHDivDofs for Faces, not implemented!");
  };
  virtual void getHDivShapes(indexType order,
                             Types::Matrix2X<prec> &shape,
                             Types::VectorX<prec> &dshape, prec xi, prec eta) {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling getHDivShapes for Faces, not implemented!");
  };

  virtual auto getHDivShapes(indexType order,
                             IntegrationPoint &IntegrationPt) -> HDivShapes {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling getHDivShapes for Faces, not implemented!");
    return {};
  };

  // L2 shapes

  virtual void getL2Dofs(std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order){};

  virtual auto getL2Shapes(indexType order,
                           IntegrationPoint &IntegrationPt) -> L2Shapes {
    return {};
  };

  virtual void getL2ShapesInternal(indexType order,
                                   Types::VectorX<prec> &shape,
                                   Types::Matrix2X<prec> &shapeDerivative,
                                   prec xsi, prec eta){};

  virtual void setL2Shapes(indexType meshId,
                           indexType order, NodeTypes type){};

  // Special Plate Shapes
  virtual void setSpecialPlateShapes(indexType meshid, indexType order,
                                     NodeTypes type) {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling setSpecialPlateShapes for Faces, not implemented!");
  };
  virtual auto getSpecialPlateDofs(indexType meshID, indexType order,
                                   NodeTypes type)
      -> std::vector<DegreeOfFreedom *> {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling setSpecialPlateShapes for Faces, not implemented!");
    return {};
  };
  virtual auto getSpecialPlateShapes(IntegrationPoint &intPoint,
                                     indexType order) -> SpecialPlateShapes {
    HierAMuS::Geometry::FacesData::throwError(
        "Error when calling getSpecialPlateShapes for Faces, not implemented!");
    return {};
  };

  virtual auto getJacobian(IntegrationPoint &IntegrationPt)
      -> Types::Matrix22<prec> = 0;

  void setAllNodeBoundaryConditionMeshId(indexType meshId,
                                         indexType dof) override;


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
  virtual void checkUpdateElement(EquationHandler &eqHandler,
                                  GeometryData &geoData) = 0;

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
  virtual auto computeMeanCoordinate()
      -> Types::Vector3<prec> = 0;

  /**
   * @brief Get the Face Normal vector
   * Only virtual method, must be implemented in the derived face classes.
   *
   * @param pointers[in], pointers to global data.
   * @return Types::Vector3<prec>, Normal vector of the face
   */
  virtual auto getFaceNormal()
      -> Types::Vector3<prec> = 0;

  virtual void set_geometry_pointers(GeometryData &geoData) = 0;

  virtual auto getVertexNumber(indexType localNumber) -> indexType = 0;
  virtual auto getVertexNumbers() -> std::vector<indexType> = 0;
  virtual void getVerts(std::vector<VertexData *> &vertsOut) = 0;

  virtual auto getEdgeNumbers() -> std::vector<indexType> = 0;
  virtual void getEdgeNumbers(std::vector<indexType> &edgesOut) = 0;
  virtual void getEdges(std::vector<EdgesData *> &edgesOut) = 0;
  virtual auto getEdgeNumber(indexType num) -> indexType = 0;

  virtual void setEdges(const std::vector<indexType> &edgesIn) = 0;

  
  virtual void setVerts(GeometryData &geoData,
                        std::vector<indexType> &vertsIn) = 0;
  virtual auto getNumberOfVerts() -> indexType = 0;

private:
  static const GeometryTypes type;
  static void throwError(const std::string &msg) {
    throw std::runtime_error(msg);
  };
};

} // namespace HierAMuS::Geometry
