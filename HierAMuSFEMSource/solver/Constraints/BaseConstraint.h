// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "Eigen/Sparse"
#include "datatypes.h"
#include "types/MatrixTypes.h"
#include <memory>

namespace HierAMuS {

class PointerCollection;
class DegreeOfFreedom;

enum class ConstraintTypes { Base = 0, GeneralLink = 1 };

struct DofRelation
{
  std::vector<DegreeOfFreedom *> m_slaveDofs;
  std::vector<DegreeOfFreedom *> m_masterDofs;
};

class BaseConstraint {
public:
  BaseConstraint();
  BaseConstraint(const BaseConstraint &other);
  virtual ~BaseConstraint() = default;

  virtual auto getCopy() -> std::shared_ptr<BaseConstraint> {
    return std::make_shared<BaseConstraint>(BaseConstraint(*this));
  };

  virtual auto getType() -> ConstraintTypes { return ConstraintTypes::Base; };

  virtual auto getDofRelation(PointerCollection &pointers) -> DofRelation {
    return DofRelation();
  };

  void setId(indexType id);
  auto getId() -> indexType;
  virtual auto getDofs(PointerCollection &pointers)
      -> std::vector<DegreeOfFreedom *> {
    return {};
  };

  virtual auto getSlaveDisplacement(PointerCollection &pointers)
      -> Types::SparseVector<prec, indexType> {
    return {};
  };

  virtual auto getSlaveNewtonSolution(PointerCollection &pointers)
      -> Types::SparseVector<prec, indexType> {
    return {};
  };
  virtual void
  modifyEquationSystem(PointerCollection &pointers,
                       Eigen::SparseMatrix<prec, 0, indexType> &Kaa,
                       Eigen::SparseMatrix<prec, 0, indexType> &Kab,
                       Eigen::SparseMatrix<prec, 0, indexType> &Kba,
                       Eigen::SparseMatrix<prec, 0, indexType> &Kbb,
                       Types::VectorX<prec> &Fa, Types::VectorX<prec> &Fb){};

  virtual void toFile(std::ofstream &out);
  virtual void fromFile(std::ifstream &in);

  virtual auto getAMatrix(PointerCollection &pointers)
      -> Types::SparseMatrix<prec, indexType> {
    return Types::SparseMatrix<prec, indexType>();};

  virtual auto getDB(PointerCollection &pointers)
      -> Types::SparseVector<prec, indexType> {
    return Types::SparseVector<prec, indexType>();
  };

private:
  indexType m_id;
};
} // namespace HierAMuS
