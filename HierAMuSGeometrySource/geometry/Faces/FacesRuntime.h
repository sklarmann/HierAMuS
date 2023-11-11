// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "geometry/Edges/EdgesData.h"
#include "geometry/VertexData.h"
#include "geometry/Faces/FaceOrientationFlags.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

#include <geometry/GeometryBaseRuntime.h>

#include <geometry/GeometryTypes.h>
#include <tuple>
#include <types/MatrixTypes.h>

#include <vector>

namespace HierAMuS {
class EquationHandler;
}

namespace HierAMuS::Geometry {
class FacesData;
class VertexRuntime;
class FacesH1Interface;
class FacesHDivInterface;

class FacesRuntime : public GeometryBaseRuntime {
public:
  FacesRuntime(FacesData &base_element);
  ~FacesRuntime() override;
  auto getGroupType() -> const GeometryTypes & override;

  virtual auto getH1Face() -> FacesH1Interface * { return NULL; };
  virtual auto getHDivFace() -> FacesHDivInterface * { return NULL; };

  virtual auto getCoordinates(IntegrationPoint &integrationPoint)
      -> Types::Vector3<prec> = 0;

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
      -> Geometry::VertexRuntime * = 0;
  virtual auto getEdge(indexType local_number)
      -> Geometry::EdgesRuntime * = 0;

  /** @brief Returns the global edge number based on the local edge number of
   * the face.
   *
   *  detailed description
   *
   *  @param[in] local_number local edge number of the face
   *  @return indexType global edge number
   */
  virtual auto getEdgeNumber(indexType local_number) -> indexType = 0;

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
    HierAMuS::Geometry::FacesRuntime::throwError(
        "Error when calling getOrientation for Faces, not implemented!");
    return {};
  };

  virtual void modifyIntegrationpoint(IntegrationPoint &IP, prec &shapeFactor,
                                      faceorientation orientation) {
    HierAMuS::Geometry::FacesRuntime::throwError(
        "Error when calling modifyInegrationpoint for Faces, not implemented!");
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
    HierAMuS::Geometry::FacesRuntime::throwError(
        "Error when calling setSpecialPlateShapes for Faces, not implemented!");
  };
  virtual auto getSpecialPlateDofs(indexType meshID, indexType order,
                                   NodeTypes type)
      -> std::vector<DegreeOfFreedom *> {
    HierAMuS::Geometry::FacesRuntime::throwError(
        "Error when calling setSpecialPlateShapes for Faces, not implemented!");
    return {};
  };
  virtual auto getSpecialPlateShapes(IntegrationPoint &intPoint,
                                     indexType order) -> SpecialPlateShapes {
    HierAMuS::Geometry::FacesRuntime::throwError(
        "Error when calling getSpecialPlateShapes for Faces, not implemented!");
    return {};
  };

  virtual auto getJacobian(IntegrationPoint &IntegrationPt) -> Types::Matrix22<prec> = 0;

  void setAllNodeBoundaryConditionMeshId(indexType meshId,
                                         indexType dof) override;

  // Paraview
  virtual void geometryToParaview(vtkPlotInterface &paraviewAdapter,
                                  indexType mainMesh, indexType subMesh){};
  virtual void computeWeightsParaview(vtkPlotInterface &paraviewAdapter,
                                      indexType mainMesh, indexType subMesh){};
  virtual void H1SolutionToParaview(vtkPlotInterface &paraviewAdapter,
                                    indexType mainMesh, indexType subMesh,
                                    indexType order,
                                    Types::VectorX<prec> &solution,
                                    std::string &name){};
  virtual void H1DataToParaview(vtkPlotInterface &paraviewAdapter,
                                indexType mainMesh, indexType subMesh,
                                Types::VectorX<prec> &Data,
                                indexType numberComponents, indexType order,
                                std::string &name){};
  virtual void projectDataToParaviewVertices(
      vtkPlotInterface &paraviewAdapter,
      indexType mainMesh, indexType subMesh, indexType order,
      IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
      indexType numberComponents, std::string name){};

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

  virtual auto getVertexNumbers() -> std::vector<indexType> = 0;
  virtual void getVertices(std::vector<VertexRuntime *> &vertsOut) = 0;

  virtual auto getEdgeNumbers() -> std::vector<indexType> = 0;
  virtual void getEdges(std::vector<EdgesRuntime *> &edgesOut) = 0;

  
  
  virtual auto getNumberOfVerts() -> indexType = 0;

    /**
   * @brief Get the global vertex number of the object using its local number.
   *
   * @param[in] localNumber, local number of the vertex
   * @return indexType, The global number of the vertex.
   */
  virtual auto getVertexNumber(indexType localNumber) -> indexType =0;

private:
  static const GeometryTypes type;
  static void throwError(const std::string &msg) {
    throw std::runtime_error(msg);
  };
};

} // namespace HierAMuS::Geometry
