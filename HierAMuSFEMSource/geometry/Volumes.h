// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <forwarddeclaration.h>
#include <geometry/Base.h>

#include <geometry/GeometryTypes.h>
#include <types/MatrixTypes.h>

#include <vector>

namespace HierAMuS::Geometry {

class Volumes : public Base {
public:
  Volumes();
  ~Volumes() override;
  auto getGroupType() -> const GeometryTypes & override;

  
  /**
   * @brief Returns the number of edges of the element.
   * @return Number of edges as an indexType
   */
  virtual auto getNumberOfEdges() -> indexType = 0;

  // Geometric mapping
  virtual void getJacobian(PointerCollection &pointers,
                           Types::Matrix33<prec> &jacobi, prec xsi, prec eta,
                           prec zeta) {
    HierAMuS::Geometry::Volumes::throwError(
        "Error when calling getJacobian for Volumes, not implemented!");
  }
  auto getJacobian(PointerCollection &pointers, IntegrationPoint &point)
      -> Types::MatrixXX<prec> override;

  // H1 Shapes
  virtual void setH1Shapes(PointerCollection &pointers, indexType meshId,
                           indexType order, NodeTypes type) {
    HierAMuS::Geometry::Volumes::throwError(
        "Error when calling setH1Shapes for Volumes, not implemented!");
  };
  virtual void setH1ShapesInternal(PointerCollection &pointers,
                                   indexType meshId, indexType order,
                                   NodeTypes type) {
    HierAMuS::Geometry::Volumes::throwError(
        "Error when calling setH1ShapesInternal for Volumes, not implemented!");
  };

  virtual void getH1Dofs(PointerCollection &pointers,
                         std::vector<DegreeOfFreedom *> &Dofs, indexType meshID,
                         indexType order) {
    HierAMuS::Geometry::Volumes::throwError(
        "Error when calling getH1Dofs for Volumes, not implemented!");
  };

  virtual void getH1DofsInternal(PointerCollection &pointers,
                                 std::vector<DegreeOfFreedom *> &Dofs,
                                 indexType meshID, indexType order) {
    HierAMuS::Geometry::Volumes::throwError(
        "Error when calling getH1DofsInternal for Volumes, not implemented!");
  };

  virtual void getH1Shapes(PointerCollection &pointers, indexType order,
                           Types::VectorX<prec> &shape,
                           Types::Matrix3X<prec> &shapeDerivative, prec xsi,
                           prec eta, prec zeta) {
    HierAMuS::Geometry::Volumes::throwError(
        "Error when calling getH1Shapes for Volumes, not implemented!");
  };

  virtual void getH1ShapesInternal(PointerCollection &pointers, indexType order,
                                   Types::VectorX<prec> &shape,
                                   Types::Matrix3X<prec> &shapeDerivative,
                                   prec xsi, prec eta, prec zeta) {
    HierAMuS::Geometry::Volumes::throwError(
        "Error when calling getH1ShapesInternal for Volumes, not implemented!");
  };

  /**
   * @brief New version of getH1Shapes.
   * Virtual method to get the H1 shapes of the element.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param order[in], order of the H1 shape functions.
   * @param IntegrationPt[in], integration point.
   * @return H1Shapes, H1 shape functions.
   */
  virtual auto getH1Shapes(PointerCollection &pointers, indexType order,
                           IntegrationPoint &IntegrationPt) -> H1Shapes {
    HierAMuS::Geometry::Volumes::throwError(
        "Function getH1Shapes for the geometry element not implemented!");
    return {};
  };
  /**
   * @brief New version of getH1ShapesInternal.
   * Computes and returns the H1 bubble functions of the volume element.
   *
   * @param pointers[in], object containing the pointers to global data.
   * @param order[in], order of the H1 shape functions.
   * @param IntegrationPt[in], integration point.
   * @return H1Shapes, H1 shape functions.
   */
  virtual auto getH1ShapesInternal(PointerCollection &pointers, indexType order,
                                   IntegrationPoint &IntegrationPt)
      -> H1Shapes {
    HierAMuS::Geometry::Volumes::throwError(
        "Function getH1Shapes for the geometry element not implemented!");
    return {};
  };

  /**
   * @brief Returns a vector of the H1 nodes for the given meshId and shape
   * function order.
   *
   * Here it is a virtual method overridden by the actual geometry object.
   *
   * @param pointers, object containtig the pointers to global data.
   * @param meshID, mesh id of the nodes for the H1 shape functions.
   * @param order, order of the H1 shape functions.
   * @return std::vector<GenericNodes *>, vector of the GenericNodes.
   */
  virtual auto getH1Nodes(PointerCollection &pointers, indexType meshID,
                          indexType order) -> std::vector<GenericNodes *> {
    return {};
  };
  virtual auto getH1NodesInternal(PointerCollection &pointers, indexType meshID,
                                  indexType order)
      -> std::vector<GenericNodes *> {
    return {};
  };

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
  virtual void checkUpdateElement(GeometryData &geoData) = 0;

private:
  static const GeometryTypes type;
  static void throwError(const std::string &msg) { throw std::runtime_error(msg); };
};

} // namespace HierAMuS
