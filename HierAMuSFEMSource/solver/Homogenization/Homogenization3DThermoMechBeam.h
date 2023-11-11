// Copyright 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "HomogeniztaionBase.h"
#include "MatrixTypes.h"
#include "control/ParameterList.h"
#include "datatypes.h"
#include "pointercollection/pointercollection.h"
#include <vector>

namespace HierAMuS {

class Homogenization3DThermoMechBeam : public HomogenizationBase {
public:
  Homogenization3DThermoMechBeam();
  Homogenization3DThermoMechBeam(const Homogenization3DThermoMechBeam &other);
  ~Homogenization3DThermoMechBeam(){};

  void init(PointerCollection &pointers, ParameterList &parameters) override;
  void computeAMatrix(PointerCollection &pointers) override;
  auto getNumberOfStrains() -> indexType override;

  auto getAMatrix() -> Types::MatrixXX<prec> & override;
  auto getDv() -> prec override;

  auto getType() -> indexType override { return 20; };

  auto getDisplacementIncrement(Types::VectorX<prec> &strainIncrement)
      -> Types::VectorX<prec> override;

  void setPeriodicDisplacements(PointerCollection& pointers,
                                Types::VectorX<prec> &strains, Types::VectorX<prec> &strainIncrement) override;


  void toFile(PointerCollection &pointers, std::ofstream &out) override;
  void fromFile(PointerCollection &pointers, std::ifstream &in) override;

private:
  void setDisplacementBoundaryConditions(PointerCollection &pointers);
  void setDisplacementBoundaryConditions2(PointerCollection &pointers);
  void setPeriodicBoundaryConditions(PointerCollection &pointers);
  void computeGeometryParameters(PointerCollection &pointers);

  
  Types::MatrixXX<prec> homogenizationMatrix;
  Types::MatrixXX<prec> homogenizationMatrixDisp;

  indexType meshIdDisp;
  indexType dispOrder;
  indexType meshIdTemp;
  indexType tempOrder;
  indexType bctype;

  std::vector<indexType> leftFacesMaster, rightFacesSlave;
  std::vector<indexType> bottomFacesMaster, topFacesSlave;
  std::vector<indexType> backFacesMaster, frontFacesSlave;

  std::vector<indexType> xEdges1, xEdges2, xEdges3, xEdges4;
  std::vector<indexType> yEdges1, yEdges2, yEdges3, yEdges4;
  std::vector<indexType> zEdges1, zEdges2, zEdges3, zEdges4;

  // corner vertices
  indexType v1, v2, v3, v4, v5, v6, v7, v8;
  

  Types::Vector3<prec> xmin, xmax;
};
} // namespace HierAMuS