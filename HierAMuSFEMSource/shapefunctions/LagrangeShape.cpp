// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <shapefunctions/LagrangeShape.h>



void LagrangeShape(prec &shape, prec &shapeDerivative, prec coor,
                   indexType order, indexType number) {

  switch (order) {
  case 1:
    switch (number) {
    case 0:
      shape = ((prec)1 - coor) / (prec)2;
      shapeDerivative = -(prec)1 / (prec)2;
      break;
    case 1:
      shape = ((prec)1 + coor) / (prec)2;
      shapeDerivative = (prec)1 / (prec)2;
      break;

    default:
      break;
    }
    break;
  case 2:
    switch (number) {
    case 0:
      shape = (coor - prec(1)) * coor / prec(2);
      shapeDerivative = (coor - prec(1)) / prec(2) + coor / prec(2);
      break;
    case 1:
      shape = (prec(1) + coor) * coor / prec(2);
      shapeDerivative = (prec(1) + coor) / prec(2) + coor / prec(2);
      break;
    case 2:
      shape = (prec(1) + coor) * (prec(1) - coor);
      shapeDerivative = (prec(1) - coor) - (prec(1) + coor);
      break;

    default:
      break;
    }
  default:
    break;
  }
}
