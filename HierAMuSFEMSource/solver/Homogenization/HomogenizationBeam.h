// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "HomogeniztaionBase.h"
#include "MatrixTypes.h"

namespace HierAMuS {

class HomogenizationBeam : public HomogenizationBase {
public:
  HomogenizationBeam();
  HomogenizationBeam(const HomogenizationBeam &other);
  ~HomogenizationBeam(){};

  void init(PointerCollection &pointers, ParameterList &parameters) override;
  void computeAMatrix(PointerCollection &pointers) override;
  auto getNumberOfStrains() -> indexType override;
      
  
  auto getAMatrix() -> Types::MatrixXX<prec> & override;
  auto getDv() -> prec override;

  auto getType() -> indexType override { return 2; };

  auto getDisplacementIncrement(Types::VectorX<prec> &strainIncrement)
      -> Types::VectorX<prec> override;

  void setPeriodicDisplacements(PointerCollection &pointers,
                                Types::VectorX<prec> &strains,
                                Types::VectorX<prec> &strainIncrement) override;



  void toFile(PointerCollection &pointers, std::ofstream &out) override;
  void fromFile(PointerCollection &pointers, std::ifstream &in) override;

private:
  void setDisplacementBoundaryConditions(PointerCollection &pointers);
  void setDisplacementBoundaryConditions2(PointerCollection &pointers);
  void setPeriodicBoundaryConditions(PointerCollection &pointers);


  Types::MatrixXX<prec> homogenizationMatrix;

  indexType meshIdDisp;
  indexType dispOrder;
  indexType bctype;

  std::vector<indexType> leftFacesMaster, rightFacesSlave;

  Types::Vector3<prec> xmin, xmax;
};
}