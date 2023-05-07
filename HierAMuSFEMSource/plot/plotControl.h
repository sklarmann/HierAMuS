// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "datatypes.h"

#include "pointercollection/pointercollection.h"

#include <memory>

namespace HierAMuS {

  class vtkPlotInterface;

  class PlotControl {
    public:
    PlotControl();
    ~PlotControl(){};

    void initialize(PointerCollection &pointers);
    void initializeMesh(PointerCollection &pointers);
    void timeUpdate(prec time);
    void toFile(PointerCollection &pointers);
    private:
    std::shared_ptr<vtkPlotInterface> m_plot;
    bool m_init, m_initMesh;
  };
}