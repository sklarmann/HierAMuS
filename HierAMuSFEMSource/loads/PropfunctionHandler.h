// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <datatypes.h>
#include <forwarddeclaration.h>

#include "Timefunction.h"
#include <Base/FEMBase.h>
#include <functional>
#include <mutex>
#include <string>
#include <vector>

namespace HierAMuS {

class PointerCollection;

class PropfunctionHandler : public FEMBase {
public:
  PropfunctionHandler();
  PropfunctionHandler(const PropfunctionHandler &other);
  ~PropfunctionHandler();
  //void setNames(PointerCollection &pointers, std::string &timename,
  //              std::string &dtname);
  //void addFunction(PointerCollection &pointers, indexType num,
  //                 std::string function, prec tmin, prec tmax);
  void addFunction(PointerCollection &pointers, indexType num,
                   std::shared_ptr<Function> &function, prec tmin, prec tmax);
  void print(PointerCollection &pointers);
  //void setDt(PointerCollection &pointers, std::string dtval);
  void set_dt(prec dtval);
  void incrTime();
  void update();
  prec getPropValue(indexType propNum);
  prec getPropValueIncr(indexType &propNum);
  prec getTime() const;
  prec getDt() const;
  void setTime(prec time);

private:
  std::vector<Timefunction> timefunctions;
  std::string timename, dtname;
  indexType numtime;
  std::mutex mutex_lock;
  prec t;
  prec dt;
};
} /* namespace HierAMuS */
