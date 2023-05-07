// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once
#include <datatypes.h>

#include <forwarddeclaration.h>
#include <types/MatrixTypes.h>
#include <vector>

namespace HierAMuS {

class Plane {
public:
  Plane() = default;
  ;
  ~Plane() = default;
  ;
  void set(const std::vector<prec> &normal, const std::vector<prec> &point);
  void set(const Types::Vector3<prec> &normal,
           const Types::Vector3<prec> &point);
  bool inPlane(const Types::Vector3<prec> &point);
  bool inPlane(prec x, prec y, prec z);

private:
  Types::Vector3<prec> normal;
  prec value;
};

} // namespace HierAMuS
