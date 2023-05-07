// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "forwarddeclaration.h"
#include <memory>

#include "macroCMD/macroCommands.h"

namespace HierAMuS {
class controlProgram {
public:
  controlProgram();
  virtual ~controlProgram()= default;;

  virtual macroCommands getMacroCommands();

  std::shared_ptr<PointerCollection> ptrCol;
};
} // namespace HierAMuS
