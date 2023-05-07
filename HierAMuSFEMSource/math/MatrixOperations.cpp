// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <types/MatrixTypes.h>
#include <datatypes.h>

#include <limits>

#include "MatrixOperations.h"

namespace HierAMuS {
namespace Math {

Types::MatrixXX<prec> AinvTimesB(Types::MatrixXX<prec> &A,
                                 Types::MatrixXX<prec> &B) {
  Types::MatrixXX<prec> retMat;
  indexType cols, rows;
  rows = A.rows();
  cols = B.cols();
  retMat.resize(rows, cols);
  auto fact = A.fullPivLu();
  // auto fact = A.jacobiSvd( Eigen::ComputeThinU | Eigen::ComputeThinV);
  // auto fact = A.bdcSvd( Eigen::ComputeThinU | Eigen::ComputeThinV);
  for (auto i = 0; i < cols; ++i) {
    retMat.block(0, i, rows, 1) = fact.solve(B.block(0, i, rows, 1));
  }
  return retMat;
}

/**
 * @brief Set the Entries To Zero Epsilon object if entry is less than
 * epsFactor*machineEpsilon.
 *
 * @tparam prec Number precision type
 * @tparam indexType Index parameter type
 * @param matrix matrix to modify
 * @param epsFactor Factor to multiply machine epsilon with, machEps*epsFactor
 * will be set to zero.
 */

void setEntriesToZeroEpsilon(Types::MatrixXX<prec> &matrix, prec epsFactor) {
  prec maxVal = 0;
  prec minVal = 0;
  for (auto i = 0; i < matrix.rows(); ++i) {
    prec max = matrix.row(i).maxCoeff();
    prec min = matrix.row(i).minCoeff();
    if (maxVal < max)
      maxVal = max;
    if (minVal > min)
      minVal = min;
  }
  if (maxVal < abs(minVal))
    maxVal = abs(minVal);
  maxVal *= std::numeric_limits<prec>::epsilon() * epsFactor;

  for (auto i = 0; i < matrix.rows(); ++i) {
    for (auto j = 0; j < matrix.cols(); ++j) {
      if (abs(matrix(i, j)) < maxVal)
        matrix(i, j) = prec(0);
    }
  }
}
} // namespace Math
} // namespace HierAMuS