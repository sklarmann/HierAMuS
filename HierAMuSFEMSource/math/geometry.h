// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once

#include <types/MatrixTypes.h>

namespace HierAMuS::Math::Geometry {

template <typename prec>
Types::Matrix33<prec> skewMatrix(const Types::Vector3<prec> &Vec) {
  Types::Matrix33<prec> retMat;
  retMat.setZero();
  retMat(0, 1) = -Vec(2);
  retMat(0, 2) = Vec(1);

  retMat(1, 0) = Vec(2);
  retMat(1, 2) = -Vec(0);

  retMat(2, 0) = -Vec(1);
  retMat(2, 1) = Vec(0);
  return retMat;
}

template <typename prec>
Types::Matrix33<prec>
getRotationMatrixBetweenTwoVectors(const Types::Vector3<prec> &VecA,
                                   const Types::Vector3<prec> &VecB) {
  auto VecAunit = VecA.normalized();
  auto VecBunit = VecB.normalized();

  prec cc = VecAunit.dot(VecBunit);
  auto rotaxis = VecAunit.cross(VecBunit);
  prec ss = rotaxis.norm();

  Types::Matrix33<prec> skewMatrix;
  skewMatrix.setZero();
  if (ss != prec(0)) {

    rotaxis.normalize();

    skewMatrix = Math::Geometry::skewMatrix(rotaxis);
  }

  Types::Matrix33<prec> result;
  result.setIdentity();

  result = result + ss * skewMatrix + (1 - cc) * skewMatrix * skewMatrix;

  return result;
}

template <typename prec>
Types::Matrix33<prec> getRotationMatrix(const Types::Vector3<prec> &Axis,
                                        prec angle) {
  auto tempAxis = Axis.normalized();
  auto skewMat = skewMatrix(tempAxis);

  Types::Matrix33<prec> result;
  result.setIdentity();

  result = result + sin(angle) * skewMat + (1 - cos(angle)) * skewMat * skewMat;

  return result;
}

template <typename prec>
Types::Matrix33<prec>
getRotationMatrix(const Types::Vector3<prec> &scaledAxis) {
  auto tempAxis = scaledAxis.normalized();
  auto skewMat = skewMatrix(tempAxis);
  prec angle = scaledAxis.norm();

  Types::Matrix33<prec> result;
  result.setIdentity();

  result = result + sin(angle) * skewMat + (1 - cos(angle)) * skewMat * skewMat;

  return result;
}

template <typename prec>
prec getRotationAngleBetweenTwoVectors(const Types::Vector3<prec> &VecA,
                                       const Types::Vector3<prec> &VecB) {
  auto VecAunit = VecA.normalized();
  auto VecBunit = VecB.normalized();

  prec cc = VecAunit.dot(VecBunit);
  prec ss = VecAunit.cross(VecBunit).norm();

  prec angle = prec(0);
  prec pi = acos(0) * prec(2);
  if (cc < prec(1) && cc > prec(-1)) {
    angle = acos(cc);
  } else if (cc >= 1) {
    angle = prec(0);
  } else if (cc <= prec(-1)) {
    angle = pi;
  }

  if (ss < prec(0)) {
    angle = prec(2) * pi - angle;
  }

  return angle;
}

} // namespace HierAMuS