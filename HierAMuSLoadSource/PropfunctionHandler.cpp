// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "PropfunctionHandler.h"

#include "Timefunction.h"


#include <mutex>

namespace HierAMuS {

PropfunctionHandler::PropfunctionHandler() : t(0), dt(0) { numtime = 0; }
PropfunctionHandler::PropfunctionHandler(const PropfunctionHandler &other)
{
  this->timefunctions = other.timefunctions;
  this->timename = other.timename;
  this->dtname = other.dtname;
  this->numtime = other.numtime;
  this->t = other.t;
  this->dt = other.dt;
};

PropfunctionHandler::~PropfunctionHandler() = default;


void PropfunctionHandler::addFunction(indexType num,
                                      std::shared_ptr<Function> &function,
                                      prec tmin, prec tmax) {
  indexType cnum = static_cast<indexType>(this->timefunctions.size());

  while (cnum < num + 1) {
    this->timefunctions.emplace_back();
    ++cnum;
  }
  this->numtime = static_cast<indexType>(cnum);
  Timefunction *temp;
  temp = &this->timefunctions[num];
  temp->set(function, tmin, tmax);
}

void PropfunctionHandler::print(spdlog::logger &Logger) {
  std::size_t num = this->timefunctions.size();
  Logger.info("Current time t:                {:>12.6e}", this->t);
  Logger.info("Current time increment dt:     {:>12.6e}\n", this->dt);
  for (auto i = 0; i < num; ++i) {
    auto temp = &this->timefunctions[i];
    
    Logger.info("Timefunction:                  {:>12}", i);
    Logger.info("Current loadfactor:            {:>12.6e}", temp->currLoad());
    Logger.info("Current loadfactor increment:  {:>12.6e}\n", temp->dcurrLoad());
  }
}


void PropfunctionHandler::set_dt(prec dtval) { this->dt = dtval; }

void PropfunctionHandler::incrTime() {
  this->t += this->dt;
  for (auto &i : this->timefunctions) {
    i.timeIncr(this->t);
  }
}

prec PropfunctionHandler::getPropValue(indexType propNum) {

  if (propNum < this->numtime) {
    this->mutex_lock.lock();
    Timefunction *temp;
    temp = &this->timefunctions[propNum];
    prec val = temp->currLoad();
    this->mutex_lock.unlock();
    return val;
  } else {
    return 0;
  }
}

prec PropfunctionHandler::getPropValueIncr(indexType &propNum) {
  if (propNum < this->numtime) {
    Timefunction *temp;
    temp = &this->timefunctions[propNum];
    return temp->dcurrLoad();
  } else {
    return 0;
  }
}

prec PropfunctionHandler::getTime() const { return this->t; }

prec PropfunctionHandler::getDt() const { return this->dt; }

void PropfunctionHandler::setTime(prec time)
{
  this->t = time;
  for (auto &i:this->timefunctions)
  {
    i.setTime(this->t);
  }
}

void PropfunctionHandler::update() {
  Timefunction *tempfun;
  std::size_t num = this->timefunctions.size();
  for (auto i = 0; i < num; ++i) {
    tempfun = &this->timefunctions[i];
    tempfun->update();
  }
}
} /* namespace HierAMuS */
