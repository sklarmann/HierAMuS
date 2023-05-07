// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once
#include "types/MatrixTypes.h"
#include <datatypes.h>

namespace HierAMuS {
void setEigenSparseMatrixZero(Types::SparseMatrix<prec, indexType> &Matrix) {
  for (int k = 0; k < Matrix.outerSize(); ++k) {
    for (Types::SparseMatrix<prec, indexType>::InnerIterator it(Matrix, k); it;
         ++it) {
      it.valueRef() = prec(0);
    }
  }
}

void setEigenSparseVectorZero(Types::SparseVector<prec, indexType> &Vector) {
  for (Types::SparseVector<prec, indexType>::InnerIterator it(Vector); it;
       ++it) {
    it.valueRef() = prec(0);
  }
}
} // namespace HierAMuS
