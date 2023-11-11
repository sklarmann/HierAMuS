// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <array>

#include <geometry/GeometryTypes.h>
#include <geometry/Volumes/VolumesDataInterface.h>

#include <types/MatrixTypes.h>

#include <vector>
#include <string>

namespace HierAMuS::Geometry {
class LinearPrismRuntime;
class LinearPrismData : public VolumesDataInterface<6,9,5,LinearPrismData> {
public:
  LinearPrismData();
  ~LinearPrismData();

  
  auto getRuntimeObject(GeometryData &geoData)
      -> std::shared_ptr<VolumesRuntime> override;

  auto getType() -> const GeometryTypes & override;

  
  
  // H1 Shapes
  void setH1Shapes(indexType meshId,
                   indexType order, NodeTypes type) override;
  void setH1ShapesInternal(indexType meshId,
                           indexType order, NodeTypes type) override;

   
  void checkUpdateElement(EquationHandler &eqHandler,
                          GeometryData &geoData) override;

  static std::string getName() { return "Linear Prism"; };

private:
  
  //indexType verts[6], edges[9], faces[5];
  static const GeometryTypes type;
};

} // namespace HierAMuS