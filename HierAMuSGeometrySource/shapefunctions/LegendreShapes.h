// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once
#include <vector>
#include <datatypes.h>



namespace HierAMuS {

struct LegendreShapeValues{
  prec shapeValue;
  prec shapeDerivative;
};

class LegendreShapes {
public:
  LegendreShapes();
  ~LegendreShapes();
  static void getShape(prec &shape, prec &shapeDeriv, prec coor, indexType functionNumber);
  static auto getShape(prec coor, indexType functionNumber) -> LegendreShapeValues; 


private:
  std::vector<std::vector<prec>> coefficients;
  std::vector<std::vector<indexType>> exponents;

};


}
