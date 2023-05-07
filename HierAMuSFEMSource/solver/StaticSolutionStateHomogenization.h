// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MatrixTypes.h"
#include "solver/Homogenization/HomogeniztaionBase.h"
#include "solver/StaticSolutionState.h"
#include <forwarddeclaration.h>



#include <memory>
#include <solver/SolutionTypes.h>

#include <Eigen/Sparse>
#include <solver/GenericSolutionState.h>
#include <vector>

#include "Homogenization/HomogeniztaionBase.h"


namespace HierAMuS {

 class PointerCollection;


class StaticSolutionStateHomogenization : public StaticSolutionState {
 public:
  StaticSolutionStateHomogenization(ParameterList &parameter);
  StaticSolutionStateHomogenization(const StaticSolutionStateHomogenization &other);
  virtual ~StaticSolutionStateHomogenization() override;

  auto getCopy() -> std::shared_ptr<GenericSolutionState> override;

  
  auto getType() -> SolutionTypes override { return SolutionTypes::StaticSolutionHomogenization; }

  void setInitialValues(indexType numberOfEquations,
                        indexType numberOfActiveEquations) override;

  
  void computeLoads(PointerCollection &pointers) override;
  void assembleSystem(PointerCollection& pointers) override;

  
  void nextSolutionStep() override;
  
  void solve(PointerCollection& pointers) override;
  auto residual() -> prec override;
  void setEquationZero() override;
  void updateSolution(PointerCollection& pointers) override;


  
  void resetSolution() override;

  
  void setStrains(PointerCollection& pointers, Types::VectorX<prec> &strains) override;
  auto getCurrentStrains() -> Types::VectorX<prec>;
  auto getStrainsIncrement() -> Types::VectorX<prec>;
  void initHomogenization(PointerCollection &pointers,
                         indexType homogenizationType,
                         ParameterList &parameters) override;
  void computeAMatrix(PointerCollection &pointers) override;
  void homogenize(PointerCollection& pointers) override;
  
  auto getHomogenizedData() -> GenericSolutionState::HomogenizedData override;

  
  void toFile(PointerCollection &pointers, std::ofstream &out) override;
  void fromFile(PointerCollection &pointers, std::ifstream &in) override;

  void RVEDatatoFile(PointerCollection &pointers, std::ofstream &out) override;
  void RVEDatafromFile(PointerCollection &pointers, std::ifstream &in) override;

  auto getCMatrix() -> Types::MatrixXX<prec>;
  auto getStresses() -> Types::VectorX<prec>;

private:
  void initHomogenizationType(PointerCollection &pointers,
                          indexType homogenizationType);
  std::unique_ptr<HomogenizationBase> homogenizationData;
  Types::VectorX<prec> currStrains, strainsIncrement;
  HomogenizedData homogenizedData;

};
} /* namespace HierAMuS */
