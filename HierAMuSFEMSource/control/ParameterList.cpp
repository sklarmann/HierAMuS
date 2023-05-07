// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "ParameterList.h"

#include <utility>
namespace HierAMuS {

prec ParameterList::getPrecVal(const std::string& name) {

  auto it = this->data.find(name);
  if(it!=this->data.end()){
    auto val = it->second(0,0);
    this->data.erase(it);
    return val;
  }

  return 0;
}
Types::MatrixXX<prec> ParameterList::getPrecMatrix(const std::string& name) {
  auto it = this->data.find(name);
  if(it!=this->data.end()){
    auto val = it->second;
    this->data.erase(it);
    return val;
  }
  return {};
}
indexType ParameterList::getIndexVal(const std::string& name) {

  auto returnVal = static_cast<indexType>(this->getPrecVal(name));
  return returnVal;
}
Types::MatrixXX<indexType> ParameterList::getIndexMatrix(const std::string& name) {

  Types::MatrixXX<prec> initMatrix = this->getPrecMatrix(name);
  indexType cols,rows;
  cols = initMatrix.cols();
  rows = initMatrix.rows();

  Types::MatrixXX<indexType> returnMatrix;
  returnMatrix.resize(rows,cols);
  for(auto i=0;i<rows;++i){
    for(auto j=0;j<cols;++j){
      returnMatrix(i,j) = static_cast<indexType>(initMatrix(i,j));
    }
  }


  return returnMatrix;
}
void ParameterList::add(const std::string &name, prec value) {
  Types::MatrixXX<prec> mat;
  mat.resize(1,1);
  mat(0,0) = value;
  this->data[name] = mat;

}
void ParameterList::add(const std::string &name, Types::MatrixXX<prec> value) {
  this->data[name] = std::move(value);

}
void ParameterList::outputRemainingValues(OutputHandler &Log) {

  for(auto & it : this->data){
    Log.all() << it.first << "=" << it.second;
  }

}
bool ParameterList::empty() {
  return this->data.empty();
}
bool ParameterList::hasParameter(std::string paramName) {
  auto it = this->data.find(paramName);
  if(it!=this->data.end()){
    return true;
  }
  return false;
}

}
