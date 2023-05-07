// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <datatypes.h>
#include <forwarddeclaration.h>

#include <functional>
#include <string>
#include <vector>
#include <memory>

namespace HierAMuS {

class Function {
public:
  void setLambda(std::function<prec(prec)> &inFunction);
  prec evaluate(prec t);

private:
	std::function<prec(prec)> function;
};

class Timefunction {
public:
  Timefunction();
  virtual ~Timefunction();
  // virtual void set(const std::string &ifunc, prec tmin, prec tmax);
  virtual void set(std::shared_ptr<Function> &function, prec tmin, prec tmax);
  virtual prec currLoad();
  virtual prec dcurrLoad();
  virtual void timeIncr(prec ctime);
  virtual void print(PointerCollection &pointers);
  virtual void update();

  virtual void setTime(prec time);

private:
  void calculate();
  bool isSet, calculated;
  std::vector<std::shared_ptr<Function>> functions;
  // std::vector<Function> functions;
  std::vector<prec> tmin, tmax;
  prec propLoadCurr;
  prec propLoadold;
  prec ctime;
};

} /* namespace HierAMuS */
