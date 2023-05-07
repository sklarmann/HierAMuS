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
#ifdef USE_MKL
#include <Eigen/PardisoSupport>
#else
#include <Eigen/SparseLU>
#include <Eigen/OrderingMethods>
#endif

namespace HierAMuS {

	
	class EigenPardisoLU : public GenericSolver {
	public:
		EigenPardisoLU() {};
		~EigenPardisoLU();
		void analyze(Eigen::SparseMatrix<prec, 0, indexType> &SpMat);
		void factorize(Eigen::SparseMatrix<prec, 0, indexType> &SpMat);
		void solve(Eigen::Matrix<prec, Eigen::Dynamic, 1> &Rhs, Eigen::Matrix<prec,
			Eigen::Dynamic, 1> &solution);
	private:
		
#ifdef USE_MKL
		Eigen::PardisoLU<Eigen::SparseMatrix<prec, 0, indexType>> solver;
#else
		Eigen::SparseLU<Eigen::SparseMatrix<prec, 0, indexType>, Eigen::AMDOrdering<indexType> > solver;
#endif // USE_MKL
		bool first;
	};


}
