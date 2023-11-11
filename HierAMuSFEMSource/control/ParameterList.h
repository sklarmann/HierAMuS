// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once
#include <map>
#include <string>
#include "datatypes.h"

#include "types/MatrixTypes.h"

namespace HierAMuS {
class OutputHandler;

class ParameterList {
public:
  ParameterList()= default;
  ~ParameterList()= default;

  void add(const std::string &name,prec value);
  void add(const std::string &name,Types::MatrixXX<prec> value);


  prec getPrecVal(const std::string& name);
  indexType getIndexVal(const std::string& name);
  Types::MatrixXX<prec> getPrecMatrix(const std::string& name);
  Types::MatrixXX<indexType> getIndexMatrix(const std::string& name);

  bool hasParameter(std::string paramName);
  bool empty();
  std::map<std::string, Types::MatrixXX<prec>>::iterator begin() {return this->data.begin();};
  std::map<std::string, Types::MatrixXX<prec>>::iterator end() {return this->data.end();};

private:
  std::map<std::string,Types::MatrixXX<prec>> data;

};

}