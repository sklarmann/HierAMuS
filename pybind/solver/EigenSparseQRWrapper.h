// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

/*
* GenericSolver.h
*
*  Created on: 02.10.2016
*      Author: Klarmann
*/

#pragma once

#include <solver/GenericSolver.h>
#include <Eigen/SparseCore>
#include<Eigen/SparseQR>
#include <Eigen/OrderingMethods>

namespace HierAMuS {

	
	class EigenSparseQR : public GenericSolver {
	public:
		EigenSparseQR() {};
		~EigenSparseQR() {};
		void analyze(Eigen::SparseMatrix<prec, 0, indexType> &SpMat);
		void factorize(Eigen::SparseMatrix<prec, 0, indexType> &SpMat);
		void solve(Eigen::Matrix<prec, Eigen::Dynamic, 1> &Rhs, Eigen::Matrix<prec,
			Eigen::Dynamic, 1> &solution);
	private:
		Eigen::SparseQR<Eigen::SparseMatrix<prec, 0, indexType>, Eigen::COLAMDOrdering<indexType> > solver;
	};

}