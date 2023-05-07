// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace HierAMuS {


	enum class SolutionTypes {
		GenericSolutionState=0,
		LinearStaticSolutionState=1,
		StaticSolutionState=2,
    StaticSolutionHomogenization = 3,
		Transient=100,
		TransientSolutionNewmark=101
	};

}