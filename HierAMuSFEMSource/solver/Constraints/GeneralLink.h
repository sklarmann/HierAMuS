// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "BaseConstraint.h"
#include "MatrixTypes.h"
#include "datatypes.h"
#include "forwarddeclaration.h"
#include "pointercollection/pointercollection.h"
#include <vector>

namespace HierAMuS {

class GeneralLink : public BaseConstraint {
public:
  GeneralLink();
  GeneralLink(const GeneralLink &other);
  ~GeneralLink() = default;

  
  auto getCopy() -> std::shared_ptr<BaseConstraint> override;

  auto getType() -> ConstraintTypes override { return ConstraintTypes::GeneralLink; };

  auto getDofs(PointerCollection &pointers)
      -> std::vector<DegreeOfFreedom *> override;

  auto getDofRelation(PointerCollection &pointers) -> DofRelation override;

  void set(PointerCollection &pointers, indexType masterDof, indexType slaveDof,
           prec factor, prec difference);

  auto getSlaveDisplacement(PointerCollection &pointers)
      -> Types::SparseVector<prec, indexType> override;
  auto getSlaveNewtonSolution(PointerCollection &pointers)
      -> Types::SparseVector<prec, indexType> override;

  void modifyEquationSystem(PointerCollection &pointers,
                            Eigen::SparseMatrix<prec, 0, indexType> &Kaa,
                            Eigen::SparseMatrix<prec, 0, indexType> &Kab,
                            Eigen::SparseMatrix<prec, 0, indexType> &Kba,
                            Eigen::SparseMatrix<prec, 0, indexType> &Kbb,
                            Types::VectorX<prec> &Fa,
                            Types::VectorX<prec> &Fb) override;

  auto getAMatrix(PointerCollection &pointers)
      -> Types::SparseMatrix<prec, indexType> override;


  void setB(prec B);

  auto getDB(PointerCollection &pointers)
      -> Types::SparseVector<prec, indexType> override;

  void toFile(std::ofstream &out) override;
  void fromFile(std::ifstream &in) override;

private:
  indexType m_masterDof, m_slaveDof;
  prec m_factor, m_difference;
};
} // namespace HierAMuS
