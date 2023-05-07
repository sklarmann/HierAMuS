//
// Created by simon on 12/16/21.
//

#include "FEMBase.h"
#include <iostream>


bool HierAMuS::FEMBase::consoleDebug = false;

namespace HierAMuS {

FEMBase::FEMBase() {
  this->consoleDebug = false;
}
void FEMBase::debug(std::string debugMessage) {
#ifndef NDEBUG
  if(this->consoleDebug){
    std::cout << "\n" << debugMessage << "\n" << std::endl;
  }
#endif
}
void FEMBase::setDebug(bool setTo) {
  this->consoleDebug = setTo;
  if(this->consoleDebug){
    std::cout << "\nEnabled Console debugging info.\n" << std::endl;
  }
}


}