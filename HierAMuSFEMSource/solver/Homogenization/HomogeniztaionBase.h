// Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MatrixTypes.h"
#include "control/ParameterList.h"
#include "datatypes.h"
#include "pointercollection/pointercollection.h"
namespace HierAMuS {

/** @brief Base class for the different homogenization schemes. It provides the
 * matrix for the relation between fixed dofs and the strains.
 *  @details The relation matrix will be a full matrix of the size of the number
 * of fixed dofs times the number of strains.
 */

class HomogenizationBase {
public:
  HomogenizationBase(){};
  HomogenizationBase(const HomogenizationBase &other){};
  virtual ~HomogenizationBase(){};

  virtual void init(PointerCollection &pointers, ParameterList &parameters){};
  virtual void computeAMatrix(PointerCollection &pointers){};
  virtual auto getAMatrix() -> Types::MatrixXX<prec> & { return dummyMatrix; };
  virtual auto getNumberOfStrains() -> indexType = 0;

  virtual auto getHomogenizationMatrix() -> Types::MatrixXX<prec> & {
    return this->dummyMatrix;
  };
  virtual auto getDv() -> prec { return prec(0.0); };

  virtual auto getType() -> indexType { return 0; };


  virtual auto getDisplacementIncrement(Types::VectorX<prec> &strainIncrement) -> Types::VectorX<prec> {
    return Types::VectorX<prec>();
  };

  virtual void
  setPeriodicDisplacements(PointerCollection& pointers,
                           Types::VectorX<prec> &strains, Types::VectorX<prec> &strainIncrement) = 0;
  
  virtual void toFile(PointerCollection &pointers, std::ofstream &out);
  virtual void fromFile(PointerCollection &pointers, std::ifstream &in);

private:
  Types::MatrixXX<prec> dummyMatrix;
};

} // namespace HierAMuS