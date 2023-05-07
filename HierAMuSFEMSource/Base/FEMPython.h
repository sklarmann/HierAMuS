//
// Created by simon on 12/19/21.
//

#pragma once
#include <forwarddeclaration.h>

namespace HierAMuS {

class FEMPython {
public:
  FEMPython();
  ~FEMPython();
  PointerCollection *getPointerCollection();
private:
  PointerCollection *ptr;
};

}