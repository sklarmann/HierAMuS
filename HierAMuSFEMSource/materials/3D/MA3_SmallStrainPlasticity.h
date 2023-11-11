// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once


#include "materials/GenericMaterialFormulation.h"
#include "types/MatrixTypes.h"

namespace HierAMuS {
namespace Materials {

class MA3_SmallStrainPlasticity : public GenericMaterialFormulation {
public:
  explicit MA3_SmallStrainPlasticity(PointerCollection *ptrCol);
  ~MA3_SmallStrainPlasticity() override;

  void readData(PointerCollection &pointers, ParameterList &list) override;
  void getMaterialData(PointerCollection &pointers,
                       MaterialTransferData &material_in_out,
                       IntegrationPoint &ip) override;

  auto getInternalVariables(PointerCollection &pointers,
                            MaterialTransferData &inoutData)
      -> std::map<std::string, Types::VectorX<prec>> override;

  auto getHistoryDataStructure(PointerCollection &pointers)
      -> const HistoryDataStructure & override;

private:
  prec m_emodul, m_nu, m_y_0, m_y_inf, m_xh, m_xd, m_eta;
  indexType m_maxIterations;
  Types::Matrix66<prec> m_material_tangent;

  const static HistoryDataStructure m_historyDataStructure;
};

} // namespace Materials
} // namespace HierAMuS