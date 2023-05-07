// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace HierAMuS {

	enum SolverTypes {
		Generic = 0,
		TypeEigenPardisoLDLT=1,
		TypeEigenPardisoLLT=2,
		TypeEigenPardisoLU=3,
		TypeEigenSimplicialLDLT=4,
		TypeEigenSimplicialLLT=5,
		TypeEigenSparseLU=6,
		TypeEigenSparseQR=7
	};

}