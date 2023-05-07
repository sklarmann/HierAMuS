// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once

#include <forwarddeclaration.h>
#include <geometry/Base.h>

#include <geometry/GeometryTypes.h>
#include <types/MatrixTypes.h>

#include <vector>

namespace HierAMuS::Geometry {

class Special : public Base {
public:
  Special();
  ~Special() override;
  auto getGroupType() -> const GeometryTypes & override;

  /**
   * @brief Returns the number of edges of the element.
   * @return Number of edges as an indexType
   */
  virtual auto getNumberOfEdges() -> indexType = 0;

  // H1 Shapes
  virtual void setH1Shapes(PointerCollection &pointers, indexType meshId,
                           indexType order, NodeTypes type) {
    this->throwError(
        "Error when calling setH1Shapes for Faces, not implemented!");
  };
  virtual void setH1ShapesInternal(PointerCollection &pointers,
                                   indexType meshId,
                                   indexType order,
                                   NodeTypes type) {
    this->throwError(
        "Error when calling setH1ShapesInternal for Faces, not implemented!");
  };

  virtual void getH1Dofs(PointerCollection &pointers,
                         std::vector<DegreeOfFreedom *> &Dofs,
                         const indexType &meshID, const indexType &order) {
    this->throwError("Error when calling getH1Dofs for Faces, not implemented!");
  };

  virtual void getH1DofsInternal(PointerCollection &pointers,
                                 std::vector<DegreeOfFreedom *> &Dofs,
                                 const indexType &meshID,
                                 const indexType &order) {
    this->throwError(
        "Error when calling getH1DofsInternal for Faces, not implemented!");
  };

  virtual void getH1Shapes(PointerCollection &pointers, const indexType &order,
                           Types::VectorX<prec> &shape,
                           Types::Matrix2X<prec> &shapeDerivative,
                           const prec &xsi, const prec &eta) {
    this->throwError(
        "Error when calling getH1Shapes for Faces, not implemented!");
  };

  virtual void getH1ShapesInternal(PointerCollection &pointers,
                                   const indexType &order,
                                   Types::VectorX<prec> &shape,
                                   Types::Matrix2X<prec> &shapeDerivative,
                                   const prec &xsi, const prec &eta) {
    this->throwError(
        "Error when calling getH1ShapesInternal for Faces, not implemented!");
  };

  virtual void setBeamVertex(indexType beamVertIn) {
    this->throwError("Error when calling setBeamVertex for Special elements, "
                     "not implemented!");
  };

private:
  static const GeometryTypes type;
  static void throwError(const std::string &msg) { throw std::runtime_error(msg); };
};

} // namespace HierAMuS
