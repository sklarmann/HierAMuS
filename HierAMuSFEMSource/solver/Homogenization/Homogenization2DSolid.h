// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "HomogeniztaionBase.h"
#include "control/ParameterList.h"
#include "pointercollection/pointercollection.h"

namespace HierAMuS {

class Homogenization2DSolid : public HomogenizationBase {
public:
  Homogenization2DSolid();
  Homogenization2DSolid(const Homogenization2DSolid &other);
  ~Homogenization2DSolid(){};

  void init(PointerCollection &pointers, ParameterList &parameters) override;
  void computeAMatrix(PointerCollection &pointers) override;
  auto getNumberOfStrains() -> indexType override;

  auto getAMatrix() -> Types::MatrixXX<prec> & override;
  auto getDv() -> prec override;

  auto getType() -> indexType override { return 1; };

  auto getDisplacementIncrement(Types::VectorX<prec> &strainIncrement)
      -> Types::VectorX<prec> override;

  void setPeriodicDisplacements(PointerCollection& pointers,
                                Types::VectorX<prec> &strains, Types::VectorX<prec> &strainIncrement) override {};
  
  void toFile(PointerCollection &pointers, std::ofstream &out) override;
  void fromFile(PointerCollection &pointers, std::ifstream &in) override;

private:
  void
  setDisplacementBoundaryConditions(PointerCollection &pointers,
                                    indexType meshIdDisp, indexType dispOrder,
                                    Types::MatrixXX<indexType> &leftedges,
                                    Types::MatrixXX<indexType> &rightedges,
                                    Types::MatrixXX<indexType> &bottomedges,
                                    Types::MatrixXX<indexType> &topedges);
  void edgeEntriesAMatrix(PointerCollection &pointers, indexType meshIdDisp,
                          indexType dispOrder,
                          Types::MatrixXX<indexType> &edges);

  void setEdgesBC(PointerCollection &pointers, indexType meshIdDisp,
                  indexType dispOrder, Types::MatrixXX<indexType> &edges);

  Types::MatrixXX<prec> homogenizationMatrix;

  indexType meshIdDisp;
  indexType dispOrder;
  indexType bctype;
  Types::MatrixXX<indexType> leftedges;
  Types::MatrixXX<indexType> rightedges;
  Types::MatrixXX<indexType> bottomedges;
  Types::MatrixXX<indexType> topedges;

  prec xmin, xmax, ymin, ymax;
};
} // namespace HierAMuS