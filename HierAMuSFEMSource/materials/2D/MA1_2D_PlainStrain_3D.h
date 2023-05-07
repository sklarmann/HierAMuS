// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <forwarddeclaration.h>

#include <materials/GenericMaterialFormulation.h>
#include <types/MatrixTypes.h>

namespace HierAMuS::Materials {
  class MA2_2D_PlainStrain_3D : public GenericMaterialFormulation {
public:
  explicit MA2_2D_PlainStrain_3D(PointerCollection *ptrCol);
  ~MA2_2D_PlainStrain_3D() override;

  void readData(PointerCollection& pointers, ParameterList &list) override;
  void getMaterialData(PointerCollection& pointers, MaterialTransferData &materialInOut, IntegrationPoint& ip) override;

    
  auto getHistoryDataStructure(PointerCollection& pointers) -> const HistoryDataStructure & override;
  
  auto getInternalVariables(PointerCollection& pointers, MaterialTransferData &inoutData) -> std::map<std::string,Types::VectorX<prec>> override;
private:
  indexType m_3D_matNum;
};

} // namespace HierAMuS