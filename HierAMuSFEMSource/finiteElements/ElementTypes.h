// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once


namespace HierAMuS::FiniteElement {
enum class Elementtypes {
  Generic,
  Edge,
  Face,
  FaceConstraint,
  Volume,
  VolumeConstraint,
  LinearPrism,
  beamInterfaceElement2D,
  beamInterfaceElement3D
};
}

