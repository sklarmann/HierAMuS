// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "materials/GenericMaterialFormulation.h"
#include "types/MatrixTypes.h"

namespace HierAMuS::Materials {
  class MA3_2D_LinearElastic_Isotrop : public GenericMaterialFormulation {
public:
  explicit MA3_2D_LinearElastic_Isotrop(PointerCollection *ptrCol);
  ~MA3_2D_LinearElastic_Isotrop() override;

  void readData(PointerCollection& pointers, ParameterList &list) override;
  void getMaterialData(PointerCollection& pointers, MaterialTransferData &materialInOut, IntegrationPoint& ip) override;

private:
  prec m_emodul{}, m_nu{}, m_thick{};
  bool m_plain_strain{};
  Types::Matrix33<prec> m_material_tangent;
};

} // namespace HierAMuS