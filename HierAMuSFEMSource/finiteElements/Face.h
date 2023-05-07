// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "GenericFiniteElement.h"
#include "MatrixTypes.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

namespace HierAMuS::FiniteElement {
class Face : public GenericFiniteElement {
  using ptrCol = PointerCollection;

public:
  Face() = default;
  ~Face() override;

  auto getType() -> Elementtypes override;

  void setFace(indexType faceIn) override;
  auto getVertexIds(PointerCollection& pointers) -> std::vector<indexType> override;

  auto getVertex(ptrCol &pointers, indexType localNumber)
      -> Geometry::Vertex & override;
  auto getEdge(ptrCol &pointers, indexType localNumber)
      -> Geometry::Edges & override;
  auto getFace(ptrCol &pointers, indexType localNumber)
      -> Geometry::Faces * override;
  

  void setAllNodeBoundaryConditionMeshId(ptrCol &pointers, indexType meshId,
                                         indexType dof) override;

  // Structure
  auto getNumberOfVertices(PointerCollection& pointers) -> indexType override;
  auto getNumberOfEdges(PointerCollection& pointers) -> indexType override;

  // Geometric mapping
  auto getJacobian(ptrCol &pointers, IntegrationPoint &IntegrationPt)
      -> Types::MatrixXX<prec> override;
  void getJacobian(ptrCol &pointers, Types::Matrix22<prec> &jacobi, prec xsi,
                   prec eta) override;

  // H1 Shapes
  void setH1Shapes(ptrCol &pointers, indexType meshid,
                   indexType order) override;
  void getH1Dofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                 indexType meshID, indexType order) override;
  auto getH1Nodes(ptrCol &pointers, indexType meshID, indexType order)
      -> std::vector<GenericNodes *> override;
  void getH1Shapes(ptrCol &pointers, indexType order,
                   Types::MatrixXX<prec> &jacobi, Types::VectorX<prec> &shape,
                   Types::MatrixXX<prec> &shapeDerivative,
                   IntegrationPoint &IntegrationPt) override;
  void getH1Shapes(ptrCol &pointers, indexType order,
                   Types::Matrix22<prec> &jacobi, Types::VectorX<prec> &shape,
                   Types::Matrix2X<prec> &dshape, prec xi, prec eta) override;

  auto getH1Shapes(ptrCol &pointers, indexType order,
                   Types::MatrixXX<prec> &jacobi,
                   IntegrationPoint &IntegrationPt)
      -> Geometry::H1Shapes override;

  auto getIntegrationPoints(ptrCol &pointers) -> IntegrationPoints override;

  // HCurl Shapes

  // Special Plate shapes
  void setSpecialPlateShapes(PointerCollection &pointers, indexType meshid,
                             indexType order) override;
  auto getSpecialPlateDofs(PointerCollection &pointers, indexType meshID,
                           indexType order)
      -> std::vector<DegreeOfFreedom *> override;
  auto getSpecialPlateShapes(PointerCollection &pointers, indexType order,
                             Types::MatrixXX<prec> &jacobi,
                             IntegrationPoint &intPoint)
      -> Geometry::SpecialPlateShapes override;

  // HDiv Shapes
  void setHDivShapes(ptrCol &pointers, indexType meshid,
                     indexType order, NodeTypes type) override;

  void getHDivDofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                   indexType meshID, indexType order) override;
  void getHDivShapes(ptrCol &pointers, indexType order,
                     Types::Matrix22<prec> &jacobi,
                     Types::Matrix2X<prec> &shape, Types::VectorX<prec> &dshape,
                     prec xi, prec eta) override;

  auto getHDivShapes(PointerCollection &pointers, indexType order,
                     Types::MatrixXX<prec> &jacobi,
                     IntegrationPoint &IntegrationPt)
      -> Geometry::HDivShapes override;

  //L2Shapes
  void getL2Dofs(ptrCol& pointers, std::vector<DegreeOfFreedom*>& Dofs,
      indexType meshID, indexType order) override {
      auto face=pointers.getGeometryData()->getFace(this->face);
      face->getL2Dofs(pointers,Dofs,meshID,order);
  };
  void setL2Shapes(ptrCol& pointers, indexType meshid,
      indexType order)override {
      auto face = pointers.getGeometryData()->getFace(this->face);
      face->setL2Shapes(pointers, meshid, order,NodeTypes::displacement);
  };
  auto getL2Shapes(ptrCol& pointers, indexType order,
      Types::MatrixXX<prec>& jacobi,
      IntegrationPoint& IntegrationPt)
      -> Geometry::L2Shapes override {
      auto face = pointers.getGeometryData()->getFace(this->face);
      return face->getL2Shapes(pointers, order,IntegrationPt);
  }
  
  void toParaviewAdapter(PointerCollection &ptrCol, vtkPlotInterface &catalyst,
                         const ParaviewSwitch &ToDo);

    
  // Paraview
  void geometryToParaview(PointerCollection &pointers,
                          vtkPlotInterface &paraviewAdapter, indexType mainMesh,
                          indexType subMesh) override;
  void computeWeightsParaview(PointerCollection &pointers,
                              vtkPlotInterface &paraviewAdapter,
                              indexType mainMesh, indexType subMesh) override;
  void H1SolutionToParaview(PointerCollection &pointers,
                            vtkPlotInterface &paraviewAdapter,
                            indexType mainMesh, indexType subMesh,
                            indexType meshId, indexType order,
                            std::string name) override;
  void projectDataToParaviewVertices(
      PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
      indexType mainMesh, indexType subMesh, indexType order,
      IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
      indexType numberComponents, std::string name) override;

protected:
  indexType face;
};
} // namespace HierAMuS::FiniteElement
