// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once

#include "materials/GenericMaterialFormulation.h"
#include "types/MatrixTypes.h"

namespace HierAMuS {
namespace Materials {

class MaterialFormulation2 : public GenericMaterialFormulation {
public:
  MaterialFormulation2(PointerCollection *ptrCol);
  ~MaterialFormulation2();
  
  void readData(PointerCollection& pointers, ParameterList &list);
  void getMaterialData(const Types::VectorX<prec> &epsilon,
                       Types::VectorX<prec> &sigma,
                       Types::MatrixXX<prec> &MaterialTangent);

private:
  prec emodul, nu;
  prec mu1, mu2, mu3;
  prec e31, e32, e33, e24, e15;
  Types::MatrixXX<prec> materialTangent;
};

} // namespace Materials
} // namespace HierAMuS