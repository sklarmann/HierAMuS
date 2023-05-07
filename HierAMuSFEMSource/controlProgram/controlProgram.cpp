// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "controlProgram.h"

#include "pointercollection/pointercollection.h"

#include "macroCMD/macroCommands.h"

namespace HierAMuS {
controlProgram::controlProgram() {
  this->ptrCol = std::make_shared<PointerCollection>();
}
	

macroCommands controlProgram::getMacroCommands() {
  return macroCommands(*this);
}

} // namespace HierAMuS
