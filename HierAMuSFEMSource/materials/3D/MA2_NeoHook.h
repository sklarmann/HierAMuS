// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <forwarddeclaration.h>

#include <materials/GenericMaterialFormulation.h>
#include <types/MatrixTypes.h>

namespace HierAMuS::Materials {

class MA2_NeoHook : public GenericMaterialFormulation {
public:
  explicit MA2_NeoHook(PointerCollection *ptrCol);
  ~MA2_NeoHook() override;

  void readData(PointerCollection& pointers, ParameterList &list) override;
  void getMaterialData(PointerCollection& pointers, MaterialTransferData &material_in_out, IntegrationPoint& ip) override;

private:
  prec m_G, m_Lambda;
};

} // namespace HierAMuS