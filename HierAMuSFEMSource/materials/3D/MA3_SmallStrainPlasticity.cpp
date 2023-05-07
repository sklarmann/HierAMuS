// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "MatrixTypes.h"
#include "materials/GenericMaterialFormulation.h"
#include <limits>

#include <materials/3D/MA3_SmallStrainPlasticity.h>

#include <control/ParameterList.h>

#include <string>

namespace HierAMuS::Materials {

MA3_SmallStrainPlasticity::MA3_SmallStrainPlasticity(PointerCollection *ptrCol)
    : GenericMaterialFormulation(ptrCol) {}

MA3_SmallStrainPlasticity::~MA3_SmallStrainPlasticity() = default;

void MA3_SmallStrainPlasticity::readData(PointerCollection& pointers, ParameterList &list) {

  this->m_emodul = list.getPrecVal("emodul");
  this->m_nu = list.getPrecVal("nu");
  this->m_y_0 = list.getPrecVal("y0");
  this->m_y_inf = list.getPrecVal("yinf");
  this->m_xh = list.getPrecVal("xh");
  this->m_xd = list.getPrecVal("xd");
  this->m_eta = list.getPrecVal("eta");
  this->m_maxIterations = list.getIndexVal("maxiterations");

  if (m_maxIterations == 0) {
    m_maxIterations = 20;
  }

  this->m_material_tangent.setZero();
  for (auto i = 0; i < 3; i++) {
    this->m_material_tangent(i, i) = prec(1) - this->m_nu;
    this->m_material_tangent(i + 3, i + 3) =
        (prec(1) - prec(2) * this->m_nu) / prec(2);
  }
  this->m_material_tangent(0, 1) = this->m_nu;
  this->m_material_tangent(1, 0) = this->m_nu;
  this->m_material_tangent(0, 2) = this->m_nu;
  this->m_material_tangent(2, 0) = this->m_nu;
  this->m_material_tangent(1, 2) = this->m_nu;
  this->m_material_tangent(2, 1) = this->m_nu;

  this->m_material_tangent *= this->m_emodul / (prec(1) + this->m_nu) /
                              (prec(1) - prec(2) * this->m_nu);
}

void MA3_SmallStrainPlasticity::getMaterialData(
  PointerCollection& pointers, MaterialTransferData &materialInOut, IntegrationPoint& ip) {

  materialInOut.materialTangent = m_material_tangent;

  // get history data
  auto ep = materialInOut.historyData->getFieldMaterialUpdateNew(0);
  auto epOld = materialInOut.historyData->getFieldMaterialUpdateOld(0);
  ep = epOld;

  auto aold = materialInOut.historyData->getFieldMaterialUpdateOld(1);
  auto anew = materialInOut.historyData->getFieldMaterialUpdateNew(1);
  prec a = aold(0, 0);
  prec a_old = a;
  // prec xE = m_emodul;
  // prec xn = m_nu;
  prec y_o = m_y_0;
  prec y_i = m_y_inf;
  prec xh = m_xh;
  prec xd = m_xd;

  Types::Matrix66<prec> cmat = m_material_tangent;

  Types::Vector6<prec> ee = materialInOut.strains - ep;
  materialInOut.stresses = m_material_tangent * ee;

  Types::Vector6<prec> sig = cmat * ee;

  prec y = y_o + xh * a + (y_i - y_o) * (prec(1) - exp(-xd * a));

  prec trsig = (sig.block(0, 0, 3, 1).sum()) / prec(3);
  Types::Vector6<prec> sd;
  sd = sig;
  sd.block(0, 0, 3, 1) -= trsig * Types::Vector3<prec>::Ones();

  prec xnt = sd(0) * sd(0) + sd(1) * sd(1) + sd(2) * sd(2) +
             prec(2) * (sd(3) * sd(3) + sd(4) * sd(4) + sd(5) * sd(5));

  prec sigv = sqrt(prec(3) * xnt / prec(2));

  prec f = sigv - y;

  if (f > y * std::numeric_limits<prec>::epsilon() * 1000) {

    prec xk = m_emodul / (prec(3) * (prec(1) - prec(2) * m_nu)); // bulk modulus
    prec xm = m_emodul / (prec(2) * (prec(1) + m_nu)); // shear modulus
    // prec eta = m_eta;
    indexType nmax = m_maxIterations;
    prec etaddt = prec(0);

    prec gam = prec(0);
    prec an = a;

    indexType nit = nmax;
    if (abs(xd) < std::numeric_limits<prec>::epsilon() * 1000)
      nit = 1;
    bool iterate = true;
    while (iterate) {
      nit--;
      prec dgam =
          f / (prec(3) * xm + xh + (y_i - y_o) * xd * exp(-xd * a) + etaddt);
      gam += dgam;
      a = an + gam;
      y = y_o + xh * a + (y_i - y_o) * (prec(1) - exp(-xd * a));
      f = sigv - (prec(3) * xm + etaddt) * gam - y;
      if (abs(f) < std::numeric_limits<prec>::epsilon() * 1000 * y || nit <= 0)
        iterate = false;
    }
    if (a < a_old)
      a = a_old;
    prec dyda = xh + (y_i - y_o) * xd * exp(-xd * a) + etaddt;
    prec beta = prec(1) - prec(3) * xm * gam / sigv;
    prec beta1 = prec(2) * xm * beta;
    prec beta2 =
        prec(2) * xm *
        (prec(1) / (prec(1) + dyda / (prec(3) * xm)) + beta - prec(1)) / xnt;

    cmat(0, 0) = xk + beta1 * prec(2) / prec(3);
    cmat(0, 1) = xk - beta1 * prec(1) / prec(3);
    cmat(0, 2) = cmat(0, 1);
    cmat(1, 0) = cmat(0, 1);
    cmat(1, 1) = cmat(0, 0);
    cmat(1, 2) = cmat(0, 1);
    cmat(2, 0) = cmat(0, 2);
    cmat(2, 1) = cmat(1, 2);
    cmat(2, 2) = cmat(0, 0);
    cmat(3, 3) = beta1 * prec(0.5);
    cmat(4, 4) = cmat(3, 3);
    cmat(5, 5) = cmat(3, 3);

    for (auto i = 0; i < 6; ++i) {
      for (auto j = 0; j <= i; ++j) {
        cmat(i, j) = cmat(i, j) - beta2 * sd(i) * sd(j);
        cmat(j, i) = cmat(i, j);
      }
    }
    prec fac = prec(3) / prec(2) * gam / sigv;
    sig.block(0, 0, 3, 1) =
        beta * sd.block(0, 0, 3, 1) + trsig * Types::Vector3<prec>::Ones();
    ep.block(0, 0, 3, 1) += fac * sd.block(0, 0, 3, 1);
    fac *= prec(2);
    sig.block(3, 0, 3, 1) = beta * sd.block(3, 0, 3, 1);
    ep.block(3, 0, 3, 1) += fac * sd.block(3, 0, 3, 1);

    trsig = (sig(0) + sig(1) + sig(2)) / prec(3);
    sd = sig;
    sd.block(0, 0, 3, 1) -= trsig * Types::Vector3<prec>::Ones();
    xnt = sd(0) * sd(0) + sd(1) * sd(1) + sd(2) * sd(2) +
          prec(2) * (sd(3) * sd(3) + sd(4) * sd(4) + sd(5) * sd(5));
    sigv = sqrt(prec(3) * xnt / prec(2));
    materialInOut.materialTangent = cmat;
    materialInOut.stresses = sig;

    anew(0, 0) = a;
  }
}


auto MA3_SmallStrainPlasticity::getInternalVariables(
  PointerCollection& pointers, MaterialTransferData &inoutData)
-> std::map<std::string, Types::VectorX<prec>> {
  std::map<std::string, Types::VectorX<prec>> internalVariables;

  auto ep = inoutData.historyData->getFieldMaterialUpdateNew(0);
  auto anew = inoutData.historyData->getFieldMaterialUpdateNew(1);
  internalVariables["PlasticStrains"] = ep;
  internalVariables["EquivalenPlasticStrain"] = anew;
  return internalVariables;
}

auto MA3_SmallStrainPlasticity::getHistoryDataStructure(PointerCollection& pointers)
-> const HistoryDataStructure & {
  return m_historyDataStructure;
}

const HistoryDataStructure
    MA3_SmallStrainPlasticity::m_historyDataStructure({}, {{6, 1}, {1, 1}});

} // namespace HierAMuS::Materials
