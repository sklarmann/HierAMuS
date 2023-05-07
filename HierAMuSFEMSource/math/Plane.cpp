// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include <math/Plane.h>

namespace HierAMuS {

void Plane::set(const std::vector<prec> &normal,
                const std::vector<prec> &point) {
  if (normal.size() != 3) {
    // TODO throw exception
  }
  if (point.size() != 3) {
    // TODO throw exception
  }
  this->value = prec(0);
  for (auto i = 0; i < 3; ++i) {
    this->normal(i) = normal[i];
    this->value += normal[i] * point[i];
  }
}

void Plane::set(const Types::Vector3<prec> &normal,
                const Types::Vector3<prec> &point)
{
  this->normal = normal;
  this->value = normal.dot(point);
}

bool Plane::inPlane(const Types::Vector3<prec> &point) {
  bool ret = false;
  if (point.size() != 3) {
    // TODO throw exception
  }

  prec val2;
  val2 = this->normal.dot(point);

  prec diff;
  diff = this->value - val2;
  prec comp;
  comp = this->value;
  if (diff < prec(0))
    diff = -diff;
  if (comp < prec(0))
    comp = -comp;
  if (diff <= comp * std::numeric_limits<prec>::epsilon()*prec(1000))
    ret = true;
  return ret;
}

bool Plane::inPlane(prec x, prec y, prec z)
{
  Types::Vector3<prec> point;
  point(0) = x;
  point(1) = y;
  point(2) = z;
  return this->inPlane(point);
}

} // namespace HierAMuS
