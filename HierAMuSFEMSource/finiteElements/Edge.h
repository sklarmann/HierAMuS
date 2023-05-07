// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <forwarddeclaration.h>


#include <finiteElements/GenericFiniteElement.h>

#include <Eigen/Dense>

#include <array>

template <class bla> class vtkSmartPointer;

 class PointerCollection;

class vtkCell;

namespace HierAMuS::FiniteElement {

class Edge : public GenericFiniteElement {
  using ptrCol = PointerCollection;

public:
  Edge()= default;;
  ~Edge() override;
  auto getElementType() -> Elementtypes override {
    return Elementtypes::Edge;
  };

  void setAllNodeBoundaryConditionMeshId(ptrCol &pointers,
                                         indexType meshId,
                                         indexType dof) override;

  auto getA1Vector(ptrCol &pointers, IntegrationPoint &integration_point)
	-> Types::Vector3<prec> override;
  
  void setEdges(std::vector<indexType> &edgesIn) override;
  auto getType() -> Elementtypes override { return Elementtypes::Edge; };
  auto getVertexId(ptrCol &pointers, indexType num) -> indexType override;
  // std::vector<indexType> getVertexIds();
  auto getVertex(ptrCol &pointers, indexType localNumber) -> Geometry::Vertex & override;
  auto getEdge(ptrCol &pointers, indexType localNumber) -> Geometry::Edges & override;
  

  // Integration Points
  auto getIntegrationPoints(ptrCol &pointers) -> IntegrationPoints override;

  // Structure
  auto getNumberOfVertices(PointerCollection& pointers) -> indexType override { return 2; };
  auto getNumberOfEdges(PointerCollection& pointers) -> indexType override { return 1; };

  // Geometric mapping
  void getJacobian(ptrCol &pointers, prec &jacobi, prec xsi) override;
  auto getJacobian(ptrCol &pointers, IntegrationPoint &IntegrationPt)
      -> Types::MatrixXX<prec> override;


  // H1 Shape Functions
  void setH1Shapes(ptrCol &pointers, indexType meshID, indexType order) override;
  void getH1Dofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                 indexType meshID, indexType order) override;
  void getH1Shapes(ptrCol &pointers, indexType order, prec jacobi,
                   Types::VectorX<prec> &shape,
                   Types::VectorX<prec> &shapeDerivative, prec xsi) override;

  auto getH1Shapes(ptrCol &pointers, indexType order,
                           Types::MatrixXX<prec> &jacobi,
                           IntegrationPoint &IntegrationPt)
      -> Geometry::H1Shapes override;

  
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

private:
  indexType m_edge;
};

} // namespace HierAMuS
