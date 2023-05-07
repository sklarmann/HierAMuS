// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <forwarddeclaration.h>

#include <finiteElements/GenericFiniteElement.h>

#include <types/MatrixTypes.h>

namespace HierAMuS {
namespace FiniteElement {

class LinearPrism : public GenericFiniteElement {

public:
  using ptrCol = PointerCollection;
  LinearPrism(){};
  ~LinearPrism();
  virtual auto getType() -> Elementtypes { return Elementtypes::LinearPrism; };

  void setVolume(indexType volIn);
  auto getVertexIds(PointerCollection& pointers) -> std::vector<indexType>;
  auto getVertex(ptrCol &pointers, indexType localNumber) -> Geometry::Vertex &;

  auto getEdge(ptrCol &pointers, indexType localNumber) -> Geometry::Edges &;
  

  // Integration points
  void getGaussPoints(indexType number, std::vector<prec> &weight,
                      std::vector<prec> &xsi, std::vector<prec> &eta,
                      std::vector<prec> &zeta);

  // Structure
  auto getNumberOfVertices(PointerCollection& pointers) -> indexType { return 6; };
  auto getNumberOfEdges(PointerCollection& pointers) -> indexType { return 9; };
  auto getNumberOfFaces(PointerCollection& pointers) -> indexType { return 5; };

  // Geometric mapping
  void getJacobian(ptrCol &pointers, Types::Matrix33<prec> &jacobi,
                   prec xsi, prec eta, prec zeta);

  // H1 Shapes
  void setH1Shapes(ptrCol &pointers, indexType meshid, indexType order);
  void getH1Dofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                 indexType meshID, indexType order);
  void getH1Shapes(ptrCol &pointers, indexType order,
                   Types::Matrix33<prec> &jacobi, Types::VectorX<prec> &shape,
                   Types::Matrix3X<prec> &dshape, prec xi,
                   prec eta, prec zeta);

  // H1 Beam Shapes
  /**
   * @brief Sets the necessary degrees of freedom based on the selected order.
   *
   * @param pointers[in] Pointers to global data.
   * @param meshid[in] Meshid for the degrees of to create.
   * @param order[in] Order of the shape functions.
   */
  void setH1BeamShapes(ptrCol &pointers, indexType meshid, indexType order);
  /**
   * @brief Gets the pointers to the degrees of freedom for the selected meshID.
   *
   * @param pointers[in] Pointers to global data.
   * @param Dofs [out] Vector with the pointers to the degrees of freedom.
   * @param meshID[in] MeshID of the degrees of freedom.
   * @param order[in] Oder of the shape functions for the selected meshID.
   */
  void getH1BeamDofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                     indexType meshID, indexType order);
  /**
   * @brief Get the beam shapes and derivatives based on the order and
   * parametric position.
   *
   * @param pointers[in] Pointers to global data.
   * @param order[in] Order of the shape functions.
   * @param jacobian[in] Value of the jacobian to map shape derivatives.
   * @param shape[out] Vector with the shape functions.
   * @param shapeDerivatives[out] Vector with the derivatives of the shape
   * functions with respect to the arc length.
   * @param xi[in] Position in parameter space.
   */
  void getH1BeamShapes(ptrCol &pointers, indexType order, prec jacobian,
                       Types::VectorX<prec> &shape,
                       Types::VectorX<prec> &shapeDerivatives, prec xi);

  /**
   * @brief Computes the jacobian of the beam element.
   *
   * @param pointers[in] Pointers to global data.
   * @return Scalar value of the jacobian.
   */
  auto getBeamJacobian(ptrCol &pointers, prec xi) -> prec;

  /**
   * @brief Returns a pointer to the edge representing the beam axis.
   *
   * @param pointers[in] Pointers to global data.
   * @return Pointer to the beam edge.
   */
  auto getBeamEdge(ptrCol &pointers) -> Geometry::Edges &;

  auto getStartTriad(ptrCol &pointers) -> Types::Matrix33<prec>;
  auto getEndTriad(ptrCol &pointers) -> Types::Matrix33<prec>;

  // Paraview
  void toParaviewAdapter(PointerCollection &ptrCol, vtkPlotInterface &catalyst,
                         const ParaviewSwitch &ToDo);
  void projectOnVertsParaview(PointerCollection &ptrCol,
                              vtkPlotInterface &catalyst,
                              Types::VectorX<prec> &values, prec &xsi,
                              prec &eta, prec &zeta, prec &weight,
                              std::string name);

private:
  indexType vol;
};

} // namespace FiniteElement
} // namespace HierAMuS
