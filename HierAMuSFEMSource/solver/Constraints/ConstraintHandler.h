// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "Eigen/Sparse"
#include "datatypes.h"
#include "types/MatrixTypes.h"
#include "BaseConstraint.h"
#include "geometry/GeometryTypes.h"
#include <vector>
#include <memory>

namespace HierAMuS {

class PointerCollection;
class DegreeOfFreedom;


class ConstraintHandler {
public:
  ConstraintHandler();
  ConstraintHandler(const ConstraintHandler &other);
  ~ConstraintHandler() = default;

  void initialize(PointerCollection &pointers);

  auto getAMatrix() -> Types::SparseMatrix<prec, indexType> &;
  auto getATranspose() -> Types::SparseMatrix<prec, indexType> &;

  auto getNewConstraint(ConstraintTypes constraintType)
      -> std::shared_ptr<BaseConstraint>;
  auto getConstraint(indexType constraintNumber)
      -> std::shared_ptr<BaseConstraint>;
  auto getNumberOfConstraints() -> indexType;

  void modifyEquationSystem(PointerCollection &pointers,
                            Eigen::SparseMatrix<prec, 0, indexType> &Kaa,
                            Eigen::SparseMatrix<prec, 0, indexType> &Kab,
                            Eigen::SparseMatrix<prec, 0, indexType> &Kba,
                            Eigen::SparseMatrix<prec, 0, indexType> &Kbb,
                            Types::VectorX<prec> &Fa,
                            Types::VectorX<prec> &Fb, bool symmetric, bool upper);

  auto getSlaveNewtonSolution(PointerCollection &pointers)
      -> Types::SparseVector<prec, indexType>;



  void GeneralLinkGeo(PointerCollection &pointers,
                      Geometry::GeometryTypes geoType,
                      const std::vector<indexType> &masterNumbers,
                      const std::vector<indexType> &slaveNumbers, indexType meshId, indexType order, indexType masterDof, indexType slaveDof, prec factor, prec difference, bool
                      reorient);



  void toFile(std::ofstream &out);
  void fromFile(std::ifstream &in);

private:
  void updateDB(PointerCollection &pointers);

  void LinkDofs(PointerCollection &pointers,
                indexType masterDof,
                indexType slaveDof, prec factor, prec difference);

  void LinkVertex(PointerCollection &pointers, indexType meshId,
                  indexType masterVert,
                  indexType slaveVert, indexType masterDof, indexType slaveDof,
                  prec factor, prec difference);

  void LinkEdges(PointerCollection &pointers, indexType meshId,
                 indexType masterEdge, indexType slaveEdge, indexType masterDof,
                 indexType slaveDof, indexType factor, indexType difference, bool reorient);

  void LinkFaces(PointerCollection &pointers, indexType meshId,
                 indexType masterFace, indexType slaveFace, indexType masterDof,
                 indexType slaveDof, indexType factor, indexType difference, bool reorient);

  Types::SparseMatrix<prec, indexType> m_A, m_ATranspose;
  Types::SparseVector<prec, indexType> m_dB, m_NewtonInc;
  std::vector<std::shared_ptr<BaseConstraint>> m_constraints;


};
} // namespace HierAMuS
