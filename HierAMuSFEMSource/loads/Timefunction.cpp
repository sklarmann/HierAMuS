// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include <loads/Timefunction.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <pointercollection/pointercollection.h>

namespace HierAMuS {

Timefunction::Timefunction() {
  this->isSet = false;
  this->calculated = false;
  this->propLoadold = static_cast<prec>(0);
  this->propLoadCurr = static_cast<prec>(0);
  this->ctime = static_cast<prec>(0);
}

Timefunction::~Timefunction() = default;

void Timefunction::calculate() {
  std::size_t num = this->functions.size();
  bool search = true;
  std::size_t i = 0;
  while (search && i < num) {
    if (this->ctime >= this->tmin[i] && this->ctime < this->tmax[i])
      search = false;
    ++i;
  }
  --i;
  if (search) {
    this->propLoadCurr = static_cast<prec>(0);
    this->calculated = true;
  } else {
    this->propLoadCurr = this->functions[i]->evaluate(this->ctime);
    //this->propLoadCurr =
    //    pointers.getUserConstants()->process(this->functions[i].func);
    this->calculated = true;
  }
}

void Timefunction::set(std::shared_ptr<Function> &function, prec tmin,
                       prec tmax) {

  this->tmin.emplace_back(tmin);
  this->tmax.emplace_back(tmax);
  this->functions.push_back(function);
}

prec Timefunction::currLoad() {
  if (!this->calculated)
    this->calculate();
  return this->propLoadCurr;
}

void Timefunction::timeIncr(prec ctime) {
  this->ctime = ctime;
  this->calculated = false;
  this->propLoadold = this->propLoadCurr;
}

prec Timefunction::dcurrLoad() {
  if (!this->calculated)
    this->calculate();
  return this->propLoadCurr - this->propLoadold;
}

void Timefunction::print(PointerCollection &pointers) {
  //std::size_t num = this->functions.size();
  //OutputHandler &Log = pointers.getInfoData()->Log;
  //Log(LogLevel::NoLog, LogLevel::NoLog)
  //    << "Current Load Level: " << std::setw(10) << this->currLoad(pointers)
  //    << " Current Load increment: " << this->dcurrLoad(pointers) << std::endl;
  //Log.output(LogLevel::NoLog, LogLevel::NoLog, "test");
  //for (auto i = 0; i < num; ++i) {
  //  Log(LogLevel::NoLog, LogLevel::NoLog)
  //      << "tmin: " << std::setw(10) << this->tmin[i]
  //      << " tmax: " << std::setw(10) << this->tmax[i]
  //      << " function: " << this->functions[i].func << std::endl;
  //}
}

void Timefunction::update() { this->propLoadold = this->propLoadCurr; }

void Timefunction::setTime(prec time)
{
  this->ctime = time;
  this->calculated = false;
  this->calculate();
  this->propLoadold = this->propLoadCurr;
}

void Function::setLambda(std::function<prec(prec)> &inFunction) {
  this->function = inFunction;
}

prec Function::evaluate(prec t) { return this->function(t); }

} /* namespace HierAMuS */
