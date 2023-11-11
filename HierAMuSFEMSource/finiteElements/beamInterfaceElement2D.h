// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once


#include <finiteElements/GenericFiniteElementInterface.h>
#include <types/MatrixTypes.h>

#include <Eigen/Dense>
#include <vector>
#include <map>

template<class bla> class vtkSmartPointer;

class vtkCell;

namespace HierAMuS::Geometry {
class EdgesData;
}

namespace HierAMuS::FiniteElement {

class beamInterfaceElement2D : public GenericFiniteElementInterface<beamInterfaceElement2D> {
  using ptrCol = PointerCollection;

public:
  beamInterfaceElement2D() = default;
  ~beamInterfaceElement2D() override;
  auto getElementType() -> Elementtypes override;;
  void set_pointers(PointerCollection &pointers) override;

  auto getGlobalEdgeNumber(PointerCollection& pointers, indexType localEdgeNumber) -> indexType;

  void setVerts(std::vector<indexType> &vertsIn) override;
  void setEdges(std::vector<indexType> &edgesIn) override;
  void setSpecial(std::vector<indexType> &specialIn) override;
  auto getType() -> Elementtypes override {
    return Elementtypes::beamInterfaceElement2D;
  };
  auto getVertexId(ptrCol &pointers, indexType num) -> indexType override;
  // std::vector<indexType> getVertexIds();
  auto getVertex(ptrCol &pointers, indexType localNumber) -> Geometry::VertexData & override;

  // Structure
  auto getNumberOfVertices(PointerCollection& pointers) -> indexType override { return 1; };
  auto getNumberOfEdges(PointerCollection& pointers) -> indexType override { return this->edges.size(); };

  void computeGeometry(ptrCol &pointers);

  // Geometric mapping
  void getJacobianInterface(ptrCol &pointers, prec &jacobi, const prec &xsi,
                            const indexType &edgeNum);

  // Special shapes
  void setShapes(ptrCol &pointers, const indexType &mesIDDisp,
                 const indexType &mesIDWarp, const indexType &mesIDRot,
                 const indexType &order);
  void computeShapes(ptrCol &pointers, const indexType &mesIDDisp,
                     const indexType &mesIDWarp, const indexType &mesIDRot,
                     const indexType &order);
  void computeShapes2(ptrCol &pointers, const indexType &mesIDDisp,
                      const indexType &mesIDRot, const indexType &order);
  void setShapes2(ptrCol &pointers, const indexType &mesIDDisp,
                  const indexType &mesIDRot, const indexType &order);

  /**
   * @brief Computes global Warping functions based on local shape functions.
   *
   * @param[in] pointers Pointer Collection containing global Data.
   * @param[in] meshid Mesh Id of the warping degrees of freedom.
   * @param[in] order Order of the local shape functions for the computation of
   * the global ones.
   */
  void computeShapesLocalWarping(ptrCol &pointers, const indexType &meshid,
                                 const indexType &order);

  /**
   * @brief Compute global Warping functions based on global polynomial.
   *
   * @param[in] pointers Pointer Collection containing global Data.
   * @param[in] meshid Mesh Id of the warping degrees of freedom.
   * @param[in] order Order of the local shape functions for the computation of
   * the global ones.
   *
   * @warning Not implemented.
   */
  void computeShapesPolynomialWarping(ptrCol &pointers, const indexType &meshid,
                                      const indexType &order);

  /**
   * @brief Get the Warping Shape functions.
   *
   * @param[in] pointers Pointer Collection containing global Data.
   * @param[out] shapeWx Shape functions for the degrees of freedom in length
   * direction.
   * @param[out] shapeWy Shape functions for the degrees of freedom in thickness
   * direction.
   * @param[out] shapeWxy Derivatives in thickness direction of the shape
   * functions for the degrees of freedom in length direction.
   * @param[out] shapeWyy Derivatives in thickness direction of the shape
   * functions for the degrees of freedom in thickness direction.
   * @param[in] edgeNum Current local edge number.
   * @param[in] eta Current local integration point coordinate for the edge.
   * @param[in] meshId Mesh Id of the associated degrees of freedom.
   */
  void getWarpingShapes(ptrCol &pointers, Types::VectorX<prec> &shapeWx,
                        Types::VectorX<prec> &shapeWy,
                        Types::VectorX<prec> &shapeWxy,
                        Types::VectorX<prec> &shapeWyy, indexType edgeNum,
                        prec eta, indexType order, indexType meshId);

  /**
   * @brief Compute the shape functions for the displacements on the surface.
   *
   * @param[in] pointers Pointer Collection containing global Data.
   * @param[in] meshID Mesh Id of the displacement degrees of freedom.
   * @param[in] order Order of the local shape functions for the computation of
   * the global ones.
   */
  void computeSurfaceDispShapes(ptrCol &pointers, const indexType &meshID,
                                const indexType &order);
  /**
   * @brief Get the Surface Displacement shape functions.
   *
   * @param[out] shapeNu Shape functions to represent the displacements.
   * @param[out] shapeNbeta Shape functions to represent the rotations.
   */
  void getSurfaceDispShapes(Types::VectorX<prec> &shapeNu,
                            Types::VectorX<prec> &shapeNbeta);

  /**
   * @brief Get the Local Surface Disp Shapes sorted into long vectors.
   *
   * @param[in] pointers Pointer Collection containing global Data.
   * @param[out] shape
   * @param[out] dshape
   * @param[in] edgeNum
   * @param[in] eta
   * @param[in] order
   * @param[in] meshId
   */
  void getLocalSurfaceDispShapesSorted(ptrCol &pointers,
                                       Types::VectorX<prec> &shape,
                                       Types::VectorX<prec> &dshape,
                                       indexType edgeNum, prec eta,
                                       indexType order, indexType meshId);
  /**
   * @brief Get Jacobian in xi direction/length direction of the element.
   *
   * @return Value of the Jacobian, length/2.
   */
  auto getJacXi() -> prec;

  /**
   * @brief Get Jacobian in eta direction/cross-section direction of the
   * element.
   *
   * @param pointers Pointers to global data.
   * @param edgeNum Current local edge number.
   * @return Value of the Jacobian, edgelength/2.
   */
  auto getJacEta(ptrCol &pointers, indexType edgeNum) -> prec;

  void setNodeMapDisp(ptrCol &pointers, const indexType &meshId);
  void setNodeMapWarp(ptrCol &pointers, const indexType &meshId);

  void setDofsOnSolid(ptrCol &pointers, const indexType &meshID,
                      const indexType &order);
  void setDofsOnVert(ptrCol &pointers, const indexType &meshID);
  void getDofsOnSolid(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                      const indexType &meshID, indexType order);
  void getDofsOnVert(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                     const indexType &meshID);
  void getNodesOnSolid(ptrCol &pointers, std::vector<GenericNodes *> &Nodes,
                       const indexType &meshID);
  void getNodesOnVert(ptrCol &pointers, std::vector<GenericNodes *> &Nodes,
                      const indexType &meshID);

  auto getNumberOfSolidNodes(ptrCol &pointers, const indexType &meshID)
      -> indexType;

  void setBCOnAllNodesSolid(ptrCol &pointers, const indexType &meshID,
                            const indexType &dof);
  void setBCOnNodeSolid(ptrCol &pointers, const indexType &meshID,
                        const indexType &node, const indexType &dof);
  void setBCOnVert(ptrCol &pointers, const indexType &meshID,
                   const indexType &dof);

  // H1 Shape Functions
  void getShapes(ptrCol &pointers, const indexType &order,
                 const indexType &meshIDDisp, const indexType &meshIDWarp,
                 const indexType &meshIDRot, Types::VectorX<prec> &shapeX,
                 Types::Matrix2X<prec> &dshapeX, Types::VectorX<prec> &shapeY,
                 Types::Matrix2X<prec> &dshapeY, prec &detj, const prec &xi,
                 const prec &eta, const indexType &localedgeNumber);
  void setH1Shapes(ptrCol &pointers, indexType meshID, indexType order) override;
  void getDofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
               const indexType &meshIDDisp, const indexType &meshIDWarp,
               const indexType &meshIDRot, const indexType &order);
  void getDofs2(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                const indexType &meshIDDisp, const indexType &meshIDRot,
                const indexType &order);
  void getH1Dofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                 indexType meshID, indexType order) override;
  void getH1Shapes(ptrCol &pointers, indexType order, prec jacobi,
                   Types::VectorX<prec> &shape,
                   Types::VectorX<prec> &shapeDerivative, prec xsi);

  auto getZCoordinate(ptrCol &pointers, const indexType &edgeNum,
                      const prec &eta) -> prec;
  auto getDA(ptrCol &pointers, const indexType &edgeNum, const prec &xi,
             const prec &eta) -> prec;
  auto getNHatBar() -> Types::Matrix2X<prec> { return this->paramsWarpingX; };
  void getXiShape(ptrCol &pointers, const indexType &order, const prec &xi,
                  Types::VectorX<prec> &shape, Types::VectorX<prec> &dshape);
  void getEtaShape(ptrCol &pointers, const indexType &meshID,
                   const indexType &edgeNum, const indexType &order,
                   const prec &eta, Types::VectorX<prec> &shape,
                   Types::VectorX<prec> &dshape);
  

  void toParaviewAdapter(ptrCol &pointers, vtkPlotInterface &catalyst,
                         const ParaviewSwitch &ToDo);
  // void getElementsLocalNodalReactions(PointerCollection& ptrCol,
  // std::map<indexType, std::vector<prec>>& vReacs);
  void setMeshIdDisp(indexType MeshIdDisp);
  void setMeshIdWarp(indexType MeshIdWarp);
  void setMeshIdRot(indexType MeshIdRot);

  void computeWarpingShapesNew(PointerCollection &pointers,
                               indexType localOrder);

  auto getInterfaceGeoElem(PointerCollection &pointers) -> Geometry::BeamInterface2D *;

  void getLocalWarpingShapesA1(PointerCollection &pointers,
                               Types::VectorX<prec> &shapes,
                               Types::VectorX<prec> &shapeDerivative,
                               indexType localEdgeNumber, prec eta);

  void getLocalWarpingShapesA2(PointerCollection &pointers,
                               Types::VectorX<prec> &shapes,
                               Types::VectorX<prec> &shapeDerivative,
                               indexType localEdgeNumber, prec eta);



  auto getA1(PointerCollection &pointers) -> Types::Vector3<prec>;
  auto getA2(PointerCollection &pointers) -> Types::Vector3<prec>;

  auto getThickness(PointerCollection &pointers) -> prec;

  auto getEdge(PointerCollection &pointers, indexType localNumber) -> Geometry::EdgesData & override;

  auto getIntegrationPoints(ptrCol &pointers) -> IntegrationPoints override;

private:
  indexType interfaceElementNumber;

  indexType vert;
  prec thickness;
  std::vector<indexType> edges;

  Types::Vector3<prec> normal;
  Types::Vector3<prec> tangent;
  Types::Vector3<prec> refCoordinates;

  std::map<indexType, indexType> dofmapWarp;
  std::map<indexType, indexType> dofmapDisp;

  std::map<indexType, indexType> NodeMapDisp;
  std::map<indexType, indexType> NodeMapWarp;

  Types::Matrix2X<prec> paramsWarpingX;
  Types::VectorX<prec> paramsWarpingY;

  Types::VectorX<prec> Nu, Nbeta;


  indexType meshIdDisp, meshIdWarp, meshIdRot;
};

} // namespace HierAMuS
