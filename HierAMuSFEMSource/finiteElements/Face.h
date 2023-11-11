// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "GenericFiniteElementInterface.h"
#include "MatrixTypes.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

namespace HierAMuS::Geometry {
class FacesData;
class FacesRuntime;
}

namespace HierAMuS::FiniteElement {
class Face : public GenericFiniteElementInterface<Face> {
  using ptrCol = PointerCollection;

public:
  Face() = default;
  ~Face() override;

  auto getType() -> Elementtypes override;
  void set_pointers(PointerCollection &pointers) override;

  void setFace(indexType faceIn) override;
  auto getVertexIds(PointerCollection& pointers) -> std::vector<indexType> override;

  auto getVertex(ptrCol &pointers, indexType localNumber)
      -> Geometry::VertexData & override;
  auto getEdge(ptrCol &pointers, indexType localNumber)
      -> Geometry::EdgesData & override;
  auto getFace(ptrCol &pointers, indexType localNumber)
      -> Geometry::FacesData * override;
  

  void setAllNodeBoundaryConditionMeshId(ptrCol &pointers, indexType meshId,
                                         indexType dof) override;

  // Structure
  auto getNumberOfVertices(PointerCollection& pointers) -> indexType override;
  auto getNumberOfEdges(PointerCollection& pointers) -> indexType override;

  // Geometric mapping
  auto getJacobian(ptrCol &pointers, IntegrationPoint &IntegrationPt)
      -> Types::Matrix22<prec>;


  // H1 Shapes
  void setH1Shapes(ptrCol &pointers, indexType meshid,
                   indexType order) override;
  void getH1Dofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                 indexType meshID, indexType order) override;
  auto getH1Nodes(ptrCol &pointers, indexType meshID, indexType order)
      -> std::vector<GenericNodes *> override;
  
  auto getH1Shapes(ptrCol &pointers, indexType order,
                   Types::Matrix22<prec> &jacobi,
                   IntegrationPoint &IntegrationPt)
      -> Geometry::H1Shapes;

  auto getIntegrationPoints(ptrCol &pointers) -> IntegrationPoints override;

  // HCurl Shapes

  // Special Plate shapes
  void setSpecialPlateShapes(PointerCollection &pointers, indexType meshid,
                             indexType order) override;
  auto getSpecialPlateDofs(PointerCollection &pointers, indexType meshID,
                           indexType order)
      -> std::vector<DegreeOfFreedom *> override;
  auto getSpecialPlateShapes(PointerCollection &pointers, indexType order,
                             Types::Matrix22<prec> &jacobi,
                             IntegrationPoint &intPoint)
      -> Geometry::SpecialPlateShapes;

  // HDiv Shapes
  void setHDivShapes(ptrCol &pointers, indexType meshid,
                     indexType order, NodeTypes type) override;

  void getHDivDofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                   indexType meshID, indexType order) override;


  auto getHDivShapes(PointerCollection &pointers, indexType order,
                     Types::Matrix22<prec> &jacobi,
                     IntegrationPoint &IntegrationPt)
      -> Geometry::HDivShapes;

  //L2Shapes
  void getL2Dofs(ptrCol& pointers, std::vector<DegreeOfFreedom*>& Dofs,
      indexType meshID, indexType order) override;;
  void setL2Shapes(ptrCol& pointers, indexType meshid,
      indexType order)override;;
  auto getL2Shapes(ptrCol &pointers, indexType order,
                   Types::Matrix22<prec> &jacobi,
                   IntegrationPoint &IntegrationPt) -> Geometry::L2Shapes;
  


    
  // Paraview
  void geometryToParaview(PointerCollection &pointers,
                          vtkPlotInterface &paraviewAdapter, indexType mainMesh,
                          indexType subMesh);
  void computeWeightsParaview(PointerCollection &pointers,
                              vtkPlotInterface &paraviewAdapter,
                              indexType mainMesh, indexType subMesh);
  void H1SolutionToParaview(PointerCollection &pointers,
                            vtkPlotInterface &paraviewAdapter,
                            indexType mainMesh, indexType subMesh,
                            indexType meshId, indexType order,
                            std::string name);
  void projectDataToParaviewVertices(
      PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
      indexType mainMesh, indexType subMesh, indexType order,
      IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
      indexType numberComponents, std::string name);

private:
  indexType m_face;
  Geometry::FacesData *m_face_object;
  std::shared_ptr<Geometry::FacesRuntime> m_face_runtime;
};
} // namespace HierAMuS::FiniteElement
