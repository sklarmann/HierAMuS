// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include <materials/MaterialFormulation2.h>

#include <control/OutputHandler.h>
#include <control/HandlingStructs.h>
#include <pointercollection/pointercollection.h>

#include <string>

#include "spdlog/fmt/ostr.h"

namespace HierAMuS {
namespace Materials {

MaterialFormulation2::MaterialFormulation2(PointerCollection *ptrCol)
    : GenericMaterialFormulation(ptrCol) {}

MaterialFormulation2::~MaterialFormulation2() {}

void MaterialFormulation2::readData(PointerCollection& pointers, ParameterList &list)
{
  this->emodul = list.getPrecVal("emodul");
  
  this->nu = list.getPrecVal("nu");
  
  this->mu1 = list.getPrecVal("mu1");
  this->mu2 = list.getPrecVal("mu2");
  this->mu3 = list.getPrecVal("mu3");

  this->e31 = list.getPrecVal("e31");
  this->e32 = list.getPrecVal("e32");
  this->e33 = list.getPrecVal("e33");
  this->e24 = list.getPrecVal("e24");
  this->e15 = list.getPrecVal("e15");

  this->materialTangent.resize(9, 9);
  this->materialTangent.setZero();

  prec tt =
      this->emodul / (prec(1) + this->nu) / (prec(1) - prec(2) * this->nu);

  this->materialTangent(0, 1) = this->nu;
  this->materialTangent(0, 2) = this->nu;
  this->materialTangent(1, 2) = this->nu;
  this->materialTangent(1, 0) = this->nu;
  this->materialTangent(2, 0) = this->nu;
  this->materialTangent(2, 1) = this->nu;
  this->materialTangent *= tt;

  this->materialTangent(0, 8) = this->e31;
  this->materialTangent(1, 8) = this->e32;
  this->materialTangent(2, 8) = this->e33;
  this->materialTangent(5, 7) = this->e24;
  this->materialTangent(4, 6) = this->e15;

  this->materialTangent(8, 0) = this->e31;
  this->materialTangent(8, 1) = this->e32;
  this->materialTangent(8, 2) = this->e33;
  this->materialTangent(7, 5) = this->e24;
  this->materialTangent(6, 4) = this->e15;

  this->materialTangent(6, 6) = this->mu1;
  this->materialTangent(7, 7) = this->mu2;
  this->materialTangent(8, 8) = this->mu3;

  for (auto i = 0; i < 3; i++) {
    this->materialTangent(i, i) = (prec(1) - this->nu) * tt;
    this->materialTangent(i + 3, i + 3) =
        ((prec(1) - prec(2) * this->nu) / prec(2)) * tt;
  }

  auto &Logger = pointers.getSPDLogger();

  Logger.info("Specified Option for Piezoelectric material:");
  Logger.info("     E-Modul:        {:>12.6e}",this->emodul);
  Logger.info("     Poisson ratio:  {:>12.6e}",this->nu);
  Logger.info("     mu1:            {:>12.6e}",this->mu1);
  Logger.info("     mu2:            {:>12.6e}",this->mu2);
  Logger.info("     mu3:            {:>12.6e}",this->mu3);
  Logger.info("     e31:            {:>12.6e}",this->e31);
  Logger.info("     e32:            {:>12.6e}",this->e32);
  Logger.info("     e33:            {:>12.6e}",this->e33);
  Logger.info("     e24:            {:>12.6e}",this->e24);
  Logger.info("     e15:            {:>12.6e}",this->e15);
  Logger.info("     Materialmatrix: \n{}",this->materialTangent);
  

}

void MaterialFormulation2::getMaterialData(
    const Types::VectorX<prec> &epsilon, Types::VectorX<prec> &sigma,
    Types::MatrixXX<prec> &MaterialTangent) {
  sigma = this->materialTangent * epsilon;
  MaterialTangent = this->materialTangent;
}

} // namespace Materials
} // namespace HierAMuS
