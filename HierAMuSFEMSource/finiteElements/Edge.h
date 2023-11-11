// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <finiteElements/GenericFiniteElementInterface.h>

#include <Eigen/Dense>

#include <array>

template <class bla> class vtkSmartPointer;

 class PointerCollection;

class vtkCell;

namespace HierAMuS::Geometry {
class EdgesRuntime;
}

namespace HierAMuS::FiniteElement {

class Edge : public GenericFiniteElementInterface<Edge> {
  using ptrCol = PointerCollection;

public:
  Edge()= default;;
  ~Edge() override;
  auto getElementType() -> Elementtypes override;

  void set_pointers(PointerCollection &pointers) override;

  void setAllNodeBoundaryConditionMeshId(ptrCol &pointers,
                                         indexType meshId,
                                         indexType dof) override;

  auto getA1Vector(ptrCol &pointers, IntegrationPoint &integration_point)
	-> Types::Vector3<prec> override;
  
  void setEdges(std::vector<indexType> &edgesIn) override;
  auto getType() -> Elementtypes override { return Elementtypes::Edge; };
  auto getVertexId(ptrCol &pointers, indexType num) -> indexType override;
  // std::vector<indexType> getVertexIds();
  auto getVertex(ptrCol &pointers, indexType localNumber) -> Geometry::VertexData & override;
  auto getEdge(ptrCol &pointers, indexType localNumber) -> Geometry::EdgesData & override;
  

  // Integration Points
  auto getIntegrationPoints(ptrCol &pointers) -> IntegrationPoints override;

  // Structure
  auto getNumberOfVertices(PointerCollection& pointers) -> indexType override { return 2; };
  auto getNumberOfEdges(PointerCollection& pointers) -> indexType override { return 1; };

  // Geometric mapping
  auto getJacobian(ptrCol &pointers, IntegrationPoint &IntegrationPt)
      -> prec;


  // H1 Shape Functions
  void setH1Shapes(ptrCol &pointers, indexType meshID, indexType order) override;
  void getH1Dofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                 indexType meshID, indexType order) override;
  
  auto getH1Shapes(ptrCol &pointers, indexType order,
                           prec jacobi,
                           IntegrationPoint &IntegrationPt)
      -> Geometry::H1Shapes;

  
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
  indexType m_edge;
  Geometry::EdgesData *m_edge_pointer;
  std::shared_ptr<Geometry::EdgesRuntime> m_edge_runtime;
};

} // namespace HierAMuS
