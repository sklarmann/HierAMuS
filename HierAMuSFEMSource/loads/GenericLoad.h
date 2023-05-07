// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <datatypes.h>



namespace HierAMuS {

	
	
	class GenericLoad {
	public:
		GenericLoad() {};
		virtual ~GenericLoad() {};
		virtual void setLoad(prec valueIn) {};
		prec getLoad() { return 0; };
	private:
	};

}