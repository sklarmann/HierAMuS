// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <forwarddeclaration.h>

#include <materials/GenericMaterialFormulation.h>
#include <types/MatrixTypes.h>

namespace HierAMuS::Materials {

class MA1_LinearElastic_Isotrop : public GenericMaterialFormulation {
public:
  explicit MA1_LinearElastic_Isotrop(PointerCollection *ptrCol);
  ~MA1_LinearElastic_Isotrop() override;

  void readData(PointerCollection& pointers, ParameterList &list) override;
  void getMaterialData(PointerCollection& pointers, MaterialTransferData &material_in_out, IntegrationPoint& ip) override;

private:
  prec m_emodul{}, m_nu{};
  Types::Matrix66<prec> m_material_tangent;
};

} // namespace HierAMuS