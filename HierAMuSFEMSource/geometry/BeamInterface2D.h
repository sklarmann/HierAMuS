// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "equations/DegreeOfFreedom.h"
#include "equations/GenericNodes.h"
#include "pointercollection/pointercollection.h"
#include <forwarddeclaration.h>
#include <geometry/Special.h>
#include <geometry/GeometryData.h>

#include <geometry/GeometryTypes.h>
#include <types/MatrixTypes.h>

#include <map>
#include <vector>

namespace HierAMuS::Geometry {

class BeamInterface2D : public Special {
public:
  BeamInterface2D();
  ~BeamInterface2D() override;
  auto getGroupType() -> const GeometryTypes & override { return HierAMuS::Geometry::BeamInterface2D::type; };

  /**
   * @brief Returns the global edge number based on the local edge number.
   *
   * @param[in] localEdgeNumber Local Edge number, zero based indexing.
   * @return indexType Global edge number.
   */
  auto getGlobalEdgeNumber(indexType localEdgeNumber) -> indexType;

  /**
   * @brief Returns the local jacobian of the edge, may be negative due to local
   * KOS
   *
   * @param[in] pointers Pointers to global data.
   * @param[in] localEdgeNumber Local edge number.
   * @param[in] eta Local parametric coordinate to evaluate at.
   * @return prec Returns the jacobian.
   */
  auto getLocalJacobian(PointerCollection &pointers, indexType localEdgeNumber,
                        prec eta) -> prec;

  auto getCoordinates(PointerCollection& pointers, IntegrationPoint& IntPoint) -> Types::Vector3<prec> override;

  /**
   * @brief Gets all degrees of freedom of the element for setup sparse matrix.
   *
   * @param[in] pointers Pointers to global data.
   * @param[out] Dofs Vector with pointers to the degrees of freedom.
   */
  void getAllDofs(PointerCollection &pointers,
                  std::vector<DegreeOfFreedom *> &Dofs);

  /**
   * @brief Sets the degrees of freedom on the surface.
   * @param[in] pointers Pointers to global data.
   * @param[in] meshId Mesh id of the degrees of freedom.
   * @param[in] order The order of the shape functions.
   * @param[in] type Type of the degrees of freedom.
   */
  void setH1ShapesSurface(PointerCollection &pointers, indexType meshId,
                          indexType order, NodeTypes type);
  /**
   * @brief Sets the degrees of freedom at the beam vertex.
   * @param[in] pointers Pointers to global data.
   * @param[in] meshId Mesh id of the degrees of freedom.
   * @param[in] type Type of the degrees of freedom.
   */
  void setDofBeamNode(PointerCollection &pointers, indexType meshId,
                      NodeTypes type);

  /**
   * @brief Returns the degrees of freedom of the solid in Dofs vector in
   * correct order to the shape functions.
   *
   * @param[in] pointers Pointers to global data.
   * @param[out] Dofs Vector of pointers to the degrees of freedom.
   * @param[in] meshID Mesh id of the degrees of freedom.
   * @param[in] localOrder Local order of the shape functions.
   */
  void getDofsOnSolid(PointerCollection &pointers,
                      std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                      indexType localOrder);
  /**
   * @brief Get all nodes at the surface with the given mesh id.
   *
   * @param[in] pointers Pointers to global data.
   * @param[out] nodes Vector of the nodes.
   * @param[in] meshId Mesh id of the degrees of freedom/nodes.
   */
  void getNodesOnSolid(PointerCollection &pointers,
                       std::vector<GenericNodes *> &nodes, indexType meshId);

  /**
   * @brief Creates a Mapping of the node ids to local numbers.
   *
   * @param[in] nodes Vector of pointers to nodes.
   * @return std::map<indexType,indexType> Map from global node ids to local
   * node ids.
   */
  auto
  getNodeMapping(std::vector<GenericNodes *> &nodes) -> std::map<indexType, indexType>;

  /**
   * @brief Returns the degrees of freedom of the beam node for the given mesh
   * id.
   *
   * @param[in] pointers Pointers to global data.
   * @param[out] Dofs Vector with degrees of freedom.
   * @param[in] meshId Mesh id of the degrees of freedom.
   */
  void getDofsOnBeamNode(PointerCollection &pointers,
                         std::vector<DegreeOfFreedom *> &Dofs,
                         indexType meshId);

  /**
   * @brief Sets the edges of the element.
   * @param[in] edgesIn Vector with the edge numbers.
   */
  void setEdges(const std::vector<indexType> &edgesIn) override;
  /**
   * @brief Sets the beam connection vertex.
   * @param[in] beamVertIn Vertex number of the connection vertex.
   */
  void setBeamVertex(indexType beamVertIn) override;
  /**
   * @brief Sets up the warping shape functions.
   * @param[in] pointers Pointers to global data.
   */
  void computeWarpingShapes(PointerCollection &pointers, indexType localOrder);
  /**
   * @brief Returns the warping shape functions for the given edge and local
   * coordinate.
   *
   * @param[in] pointers Pointers to global data.
   * @param[out] shapes VectorX of the shape functions values.
   * @param[out] shapeDerivative VectorX of the derivatives of the shape
   * functions with respect to thickness direction.
   * @param[in] localEdgeNumber The local edge number of the element.
   * @param[in] eta Local parametric coordinate to evaluation the shape
   * functions at.
   */
  void getLocalWarpingShapesA1(PointerCollection &pointers,
                               Types::VectorX<prec> &shapes,
                               Types::VectorX<prec> &shapeDerivative,
                               indexType localEdgeNumber, prec eta);

  void getLocalWarpingShapesA2(
    PointerCollection &pointers, Types::VectorX<prec> &shapes,
    Types::VectorX<prec> &shapeDerivative, indexType localEdgeNumber,
    prec eta);
  /**
   * @brief Returns to local H1 shapes and derivatives with respect to A2
   * direction.
   *
   * @param[in] pointers Pointers to global data.
   * @param[out] shapes Shape functions ordered according to the selected edge.
   * @param[out] shapeDerivative Shape function derivatives ordered according to
   * the selected edge.
   * @param[in] localEdgeNumber Local edge number to compute shapes for.
   * @param[in] eta Local parametric coordinate to evaluate shape function at.
   */
  void getLocalH1ShapesA2(PointerCollection &pointers,
                          Types::VectorX<prec> &shapes,
                          Types::VectorX<prec> &shapeDerivative,
                          indexType localEdgeNumber, prec eta,
                          indexType shapeOrder, indexType meshId);

  void getH1ShapesA1(PointerCollection &pointers, Types::VectorX<prec> &shapes,
                     Types::VectorX<prec> &shapeDerivative, prec xi,
                     indexType order);
  /**
   * @brief Transforms global edge numbers to the local ones of the element.
   * @param[in,out] edgeNumbersInOut Vector of global edge numbers as input,
   * gets transformed to local edge numbers.
   */
  void globalToLocalEdgeNumbers(std::vector<indexType> &edgeNumbersInOut);
  /**
   * @brief Returns the cross section coordinate inside the specific edge.
   * @param[in] pointers Pointers to global data.
   * @param[in] eta Parametric coordinate inside edge.
   * @param[in] edgenum Local edge number inside the element.
   * @return Cross section coordinate.
   */
  auto getCrossSectionPosition(PointerCollection &pointers, prec eta,
                               indexType edgenum) -> prec;

  /**
   *
   * @param[in] pointers Pointers to global data.
   * @param[in] localEdgeNumber Local edge number.
   * @param[in] xi Local parametric coordinate in length direction.
   * @param[in] eta Local parametric coordinate in thickness direction.
   * @return Jacobi determinant.
   */
  auto getDA(PointerCollection &pointers, indexType localEdgeNumber, prec xi,
             prec eta) -> prec;

  auto getA1() -> Types::Vector3<prec> { return this->A1; };
  auto getA2() -> Types::Vector3<prec> { return this->A2; };

  void setBCOnVert(PointerCollection &pointers, const indexType &meshID,
                   const indexType &dof);

  void setBCOnSolid(PointerCollection &pointers, indexType meshID,
                    indexType dof);

  auto getNumberOfEdges() -> indexType override { return this->edges.size(); };

  auto getThickness() -> prec { return this->thickness; };

  auto getVertex(PointerCollection &pointers) -> Vertex & {
    return pointers.getGeometryData()->getVertex(this->beamNode);
  }

  void print(PointerCollection &pointers) override {};

private:
  /**
   * @brief Computes geometric parameters of the element.
   * @param[in] pointers Pointers to global data.
   */
  void computeGeometry(PointerCollection &pointers);
  /**
   * @brief Sets up a global to local Vertex id mapping to ensure continuity of
   * warping shape functions.
   * @param[in] pointers Pointers to global data.
   */
  void setUpGlobalLocalVertexMapping(PointerCollection &pointers);

  /**
   * @brief Computes the alpha and beta parameters for the warping shape
   * functions.
   * @param[in] pointers Pointers to global data.
   */
  void computeAlphaBetaParameters(PointerCollection &pointers);

  static const GeometryTypes type;
  static void throwError(const std::string &msg) { throw std::runtime_error(msg); };

  prec thickness;
  indexType beamNode;
  indexType warpingOrder;
  bool initialized;
  Types::Vector3<prec> surfaceCoordinat, A1, A2;

  std::vector<indexType> edges;
  Types::Matrix2X<prec> alphaBeta;
  Types::VectorXT<prec> alphaA2;

  std::map<indexType, indexType> globalLocalVertIdMap;
  std::map<indexType, indexType> globalLocalEdgeIdMap;
  indexType totalShapes;
};

} // namespace HierAMuS
