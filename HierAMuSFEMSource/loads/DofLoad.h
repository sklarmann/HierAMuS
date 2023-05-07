// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <datatypes.h>

namespace HierAMuS {

	
	class DofLoad {
	public:
		DofLoad();
		~DofLoad();
		void setLoad(const prec &valueIn, const indexType &dofNumber);
		//void getLoad(prec valueOut, indexType dofNumber);
	private:
		prec value;
		indexType dof;
	};

	
	DofLoad::DofLoad() {
		this->value = 0;
	}


	
	DofLoad::~DofLoad() {

	}

	
	void DofLoad::setLoad(const prec &valueIn, const indexType &dofNumber) {
		this->value = valueIn;
		this->dof = dofNumber;
	}

	//
	//DofLoad::getLoad()
}
