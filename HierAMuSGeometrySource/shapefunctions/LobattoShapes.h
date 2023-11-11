// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once
#include <vector>
#include <datatypes.h>



namespace HierAMuS {

struct LobattoShapesValues{
  prec shapeValue;
  prec shapeDerivative;
};

class LobattoShapes {
public:
  LobattoShapes();
  ~LobattoShapes();
  static void getShape(prec &shape, prec &shapeDeriv, prec coor, indexType functionNumber);
  static auto getShape( prec coor, indexType functionNumber) -> LobattoShapesValues;


private:
  std::vector<std::vector<prec>> coefficients;
  std::vector<std::vector<indexType>> exponents;

};


}
