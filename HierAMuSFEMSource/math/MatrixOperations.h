// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "types/MatrixTypes.h"
#include "datatypes.h"


namespace HierAMuS {
namespace Math {

Types::MatrixXX<prec> AinvTimesB(Types::MatrixXX<prec> &A,
                                 Types::MatrixXX<prec> &B);

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

void setEntriesToZeroEpsilon(Types::MatrixXX<prec> &matrix, prec epsFactor);
} // namespace Math
} // namespace HierAMuS