// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "SolverTypesWrapper.h"


void HierAMuS::SolverTypesWrapper::registerFunctions()
{
  this->temp
      .value("TypeEigenPardisoLDLT",
             HierAMuS::SolverTypes::TypeEigenPardisoLDLT)
      .value("TypeEigenPardisoLLT",
             HierAMuS::SolverTypes::TypeEigenPardisoLLT)
      .value("TypeEigenPardisoLU", HierAMuS::SolverTypes::TypeEigenPardisoLU)
      .value("TypeEigenSimplicialLDLT",
             HierAMuS::SolverTypes::TypeEigenSimplicialLDLT)
      .value("TypeEigenSimplicialLLT",
             HierAMuS::SolverTypes::TypeEigenSimplicialLLT)
      .value("TypeEigenSparseLU", HierAMuS::SolverTypes::TypeEigenSparseLU)
      .value("TypeEigenSparseQR", HierAMuS::SolverTypes::TypeEigenSparseQR)

      ;

}


