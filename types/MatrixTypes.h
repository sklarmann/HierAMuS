// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#ifdef __CUDACC__
#undef __CUDACC__
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "datatypes.h"

namespace HierAMuS {
//     template<typename prec>
//     struct Matrix {
//             typename Eigen::Matrix<prec, Eigen::Dynamic, 1> VectorX;
//         };
//     template<typename tt>
namespace Types {
template <typename prec> using Vector2 = typename Eigen::Matrix<prec, 2, 1>;
template <typename prec> using Vector3 = typename Eigen::Matrix<prec, 3, 1>;
template <typename prec> using Vector6 = typename Eigen::Matrix<prec, 6, 1>;
template <typename prec> using Vector2T = typename Eigen::Matrix<prec, 1, 2>;
template <typename prec> using Vector3T = typename Eigen::Matrix<prec, 1, 3>;
template <typename prec> using Vector6T = typename Eigen::Matrix<prec, 1, 6>;

template <typename prec>
using VectorX = typename Eigen::Matrix<prec, Eigen::Dynamic, 1>;
template <typename prec>
using VectorXT = typename Eigen::Matrix<prec, 1, Eigen::Dynamic>;

template <typename prec> using Matrix22 = typename Eigen::Matrix<prec, 2, 2>;
template <typename prec> using Matrix33 = typename Eigen::Matrix<prec, 3, 3>;
template <typename prec> using Matrix66 = typename Eigen::Matrix<prec, 6, 6>;

template <typename prec>
using Matrix2X = typename Eigen::Matrix<prec, 2, Eigen::Dynamic>;
template <typename prec>
using Matrix3X = typename Eigen::Matrix<prec, 3, Eigen::Dynamic>;
template <typename prec>
using Matrix6X = typename Eigen::Matrix<prec, 6, Eigen::Dynamic>;

template <typename prec>
using MatrixXX = typename Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic>;
template <typename prec, int sx, int sy>
using Matrix = typename Eigen::Matrix<prec, sx, sy>;

/*template <typename prec>
using typedef Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> MaxtrixXX;*/
template<typename prec, typename indexType>
using SparseMatrix = typename Eigen::SparseMatrix<prec,0,indexType>;
template<typename prec, typename indexType>
using SparseVector = typename Eigen::SparseVector<prec,0,indexType>;



} // namespace Types
} // namespace HierAMuS
