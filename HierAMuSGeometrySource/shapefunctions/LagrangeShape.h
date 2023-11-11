// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once
#include <datatypes.h>

void LagrangeShape(prec &shape, prec &shapeDerivative, prec coor,
                   indexType order, indexType number);
