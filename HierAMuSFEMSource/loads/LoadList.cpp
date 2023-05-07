// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <loads/LoadList.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <loads/PropfunctionHandler.h>
#include <pointercollection/pointercollection.h>

namespace HierAMuS {

loadList::loadList() {
}

loadList::~loadList() = default;

void loadList::setLoad(indexType propNum, indexType Dof, prec loadValue,
                       bool add) {
  if (add) {
    if (this->loads.find(propNum) != this->loads.end()) {
      if (this->loads[propNum].find(Dof) != this->loads[propNum].end()) {
        this->loads[propNum][Dof] += loadValue;
      } else {
        this->loads[propNum][Dof] = loadValue;
      }
    } else {
      this->loads[propNum][Dof] = loadValue;
    }
  } else {
    this->loads[propNum][Dof] = loadValue;
  }
}

prec loadList::getLoad(PointerCollection& pointers, indexType Dof) {
  typename std::map<indexType, std::map<indexType, prec>>::iterator it;
  typename std::map<indexType, prec>::iterator it2;
  indexType propnum;
  it = this->loads.begin();
  while (it != this->loads.end()) {
    propnum = it->first;
    it2 = it->second.find(Dof);
    if (it2 != it->second.end()) {
      return it2->second *
             pointers.getPropLoads()->getPropValue(propnum);
    }
    ++it;
  }
  return 0;
}

prec loadList::getLoadIncr(PointerCollection& pointers, indexType Dof) {
  typename std::map<indexType, std::map<indexType, prec>>::iterator it;
  typename std::map<indexType, prec>::iterator it2;
  indexType propnum;
  it = this->loads.begin();
  while (it != this->loads.end()) {
    propnum = it->first;
    it2 = it->second.find(Dof);
    if (it2 != it->second.end()) {
      return it2->second *
             pointers.getPropLoads()->getPropValueIncr(propnum);
    }
    ++it;
  }
  return 0;
}

void loadList::computeLoads(PointerCollection &pointers,
                            std::vector<indexType> &ids,
                            std::vector<prec> &loads,
                            std::vector<prec> &loadincs) {

  ids.clear();
  loads.clear();
  loadincs.clear();
  typedef typename std::map<indexType, std::map<indexType, prec>>::iterator
      MapIterator1;
  MapIterator1 it1;
  typedef typename std::map<indexType, prec>::iterator MapIterator2;
  MapIterator2 it2;
  it1 = this->loads.begin();
  while (it1 != this->loads.end()) {
    indexType propnum = it1->first;
    auto solState = pointers.getSolutionState();
    auto propHandler = solState->getProps();
    prec prop = propHandler->getPropValue(propnum);
    prec dprop = propHandler->getPropValueIncr(propnum);
    it2 = it1->second.begin();
    while (it2 != it1->second.end()) {
      ids.push_back(it2->first);
      loads.push_back(prop * it2->second);
      loadincs.push_back(dprop * it2->second);
      // loads[it2->first] += prop*it2->second;
      // loadincs[it2->first] += dprop*it2->second;
      ++it2;
    }
    ++it1;
  }
}

void loadList::print(PointerCollection& pointers) {
  auto &Logger = pointers.getSPDLogger();

  Logger.info("Load list information:");

  for (auto & load : this->loads) {
    Logger.info("Loads assciated with prop function:  {:>10}",load.first);
    for (auto jt = load.second.begin(); jt != load.second.end(); ++jt) {
      Logger.info("  Load on Dof {:>12} with value: {:>12.6e}",jt->first,jt->second);
    }
  }
}

} // namespace HierAMuS
