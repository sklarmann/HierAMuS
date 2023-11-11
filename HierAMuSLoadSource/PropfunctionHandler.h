// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "datatypes.h"

#include "Timefunction.h"
#include <mutex>
#include <string>
#include <vector>


#include "spdlog/spdlog.h"

namespace HierAMuS {

class PropfunctionHandler {
public:
  PropfunctionHandler();
  PropfunctionHandler(const PropfunctionHandler &other);
  ~PropfunctionHandler();
  void addFunction(indexType num,
                   std::shared_ptr<Function> &function, prec tmin, prec tmax);
  void print(spdlog::logger &Logger);

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
