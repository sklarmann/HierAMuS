// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <forwarddeclaration.h>

#include <materials/GenericMaterialFormulation.h>
#include <types/MatrixTypes.h>

#include <memory>
#ifdef OPENMP
#include <omp.h>
#endif



namespace HierAMuS::Materials {

class MAS1_Homogenization : public GenericMaterialFormulation {
public:
  explicit MAS1_Homogenization(PointerCollection *ptrCol);
  ~MAS1_Homogenization() override;

  void readData(PointerCollection& pointers, ParameterList &list) override;
  void getMaterialData(PointerCollection& pointers, MaterialTransferData &material_in_out, IntegrationPoint& ip) override;

  void initRVE(PointerCollection &pointers, IntegrationPoint &ip) override;
  void setRVE(PointerCollection &pointers, PointerCollection &RVE) override;
  void updateRVEHistory(PointerCollection &pointers, IntegrationPoint& ip) override;

private:
  auto getFileName(std::shared_ptr<PointerCollection> pointers,
                   IntegrationPoint &ip)
      -> std::string;

  PointerCollection *m_RVE;
  indexType m_maxIterations;
  prec m_refNorm, m_refENorm;

#ifdef OPENMP
  omp_lock_t m_normLock; // OPENMP lock for norm update
#endif
};

} // namespace HierAMuS