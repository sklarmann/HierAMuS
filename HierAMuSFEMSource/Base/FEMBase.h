//
// Created by simon on 12/16/21.
//

#pragma once
#include <string>

namespace HierAMuS {

class FEMBase {
public:
  FEMBase();
  void debug(std::string debugMessage);
  void setDebug(bool setTo);


private:
  static bool consoleDebug;
};



}