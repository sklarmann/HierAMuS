// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "forwarddeclaration.h"

namespace HierAMuS {
class macroCommands {
public:
  macroCommands() = delete;
  macroCommands(controlProgram &prog);
  virtual ~macroCommands();

protected:
  controlProgram *program;
};
} // namespace HierAMuS
