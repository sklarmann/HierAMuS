// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once


#include <types/MatrixTypes.h>
#include "datatypes.h"
#include <iostream>

namespace HierAMuS {
class NodeSet;
class vtkPlotInterface;
class GenericNodes;

namespace Geometry {

struct H1Shapes {
  Types::VectorX<prec> shapes;
  Types::MatrixXX<prec> shapeDeriv;

  H1Shapes(indexType numShapes, indexType numDerivs)
      : shapes(numShapes), shapeDeriv(numDerivs, numShapes){};
  H1Shapes(H1Shapes &other)
      : shapes(other.shapes),
        shapeDeriv(other.shapeDeriv) {
            // std::cout << "copy constructor called" << std::endl;
        };
  H1Shapes(H1Shapes &&other)
      : shapes(std::move(other.shapes)),
        shapeDeriv(std::move(other.shapeDeriv)){};

  void operator=(H1Shapes &&other) {
    shapes = std::move(other.shapes);
    shapeDeriv = std::move(other.shapeDeriv);
    std::cout << "move" << std::endl;
  };
  void operator=(H1Shapes &other) {
    shapes = other.shapes;
    shapeDeriv = other.shapeDeriv;
    std::cout << "assignment called" << std::endl;
  };
};

struct L2Shapes {
  Types::VectorX<prec> shapes;
};

struct HDivShapes {
  Types::MatrixXX<prec> shapes;
  Types::MatrixXX<prec> shapeDeriv;
};

struct SpecialPlateShapes {
  Types::Matrix2X<prec> shapes;
  Types::Matrix2X<prec> shapeDx;
  Types::Matrix2X<prec> shapeDy;
};
} // namespace Geometry
} // namespace HierAMuS