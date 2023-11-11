// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <geometry/GeometryBaseRuntime.h>

#include "geometry/GeometryShape.h"
#include <geometry/GeometryTypes.h>
#include <types/MatrixTypes.h>

#include <vector>


namespace HierAMuS::Geometry {
class GeometryData;
class VolumesData;
class VertexRuntime;
class EdgesRuntime;
class FacesRuntime;
class VolumesH1Interface;

class VolumesRuntime : public GeometryBaseRuntime {
public:
  VolumesRuntime(VolumesData &data_element);
  ~VolumesRuntime() override;
  auto getGroupType() -> const GeometryTypes & override;

  virtual auto getH1Volume() -> VolumesH1Interface * { return NULL; };

  /**
   * @brief Returns the number of edges of the element.
   * @return Number of edges as an indexType
   */
  virtual auto getNumberOfEdges() -> indexType = 0;
  /**
   * @brief Get the Number Of Faces object
   *
   * @return indexType Number of faces as an indexType
   */
  virtual auto getNumberOfFaces() -> indexType = 0;


  virtual auto getCoordinates(IntegrationPoint &IntPoint)
      -> Types::Vector3<prec> = 0;

  virtual auto getJacobian(IntegrationPoint &point)
      -> Types::Matrix33<prec> = 0;


  // Paraview
  virtual void geometryToParaview(vtkPlotInterface &paraviewAdapter,
                                  indexType mainMesh, indexType subMesh){};
  virtual void computeWeightsParaview(vtkPlotInterface &paraviewAdapter,
                                      indexType mainMesh, indexType subMesh){};
  virtual void H1SolutionToParaview(vtkPlotInterface &paraviewAdapter,
                                    indexType mainMesh, indexType subMesh,
                                    indexType order, Types::VectorX<prec> &solution,
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


  virtual auto getVertexNumbers() -> std::vector<indexType> = 0;
  virtual void getVertices(std::vector<VertexRuntime *> &vertsOut) = 0;

  virtual void getEdgeNumbers(std::vector<indexType> &edgesOut) = 0;
  virtual void getEdges(std::vector<EdgesRuntime *> &edgesOut) = 0;


  virtual void getFaceNumbers(std::vector<indexType> &facesOut) = 0;

  virtual void getFaces(std::vector<FacesRuntime *> &faces) = 0;

 
  virtual auto getNumberOfVerts() -> indexType = 0;

private:
  static const GeometryTypes type;
  static void throwError(const std::string &msg) {
    throw std::runtime_error(msg);
  };
};

} // namespace HierAMuS::Geometry
