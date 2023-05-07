//
// Created by simon on 12/19/21.
//

#include "FEMPython.h"
#include <pointercollection/pointercollection.h>

namespace HierAMuS {

FEMPython::FEMPython() {
  this->ptr = new PointerCollection;
}

FEMPython::~FEMPython() {
  delete this->ptr;
}

PointerCollection *FEMPython::getPointerCollection() { return this->ptr; }
}