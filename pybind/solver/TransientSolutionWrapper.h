// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

/*
* TransientSolutionNewmark.h
*
*  Created on: 08.08.2016
*      Author: klari
*/

#pragma once

#include <forwarddeclaration.h>

#include <solver/GenericSolutionState.h>
#include <map>

#include <vector>

namespace HierAMuS {

	
	class PointerCollection;

	
	class TransientSolution : public GenericSolutionState {
	public:
		TransientSolution(PointerCollection *pointers, ParameterList &parameter) : GenericSolutionState(pointers, parameter) {};
		~TransientSolution() {};
		virtual SolutionTypes getType() { return SolutionTypes::Transient; }
    

	private:
	};
} /* namespace HierAMuS */
