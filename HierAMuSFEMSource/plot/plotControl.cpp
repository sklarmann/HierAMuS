// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "plotControl.h"
#include "control/HandlingStructs.h"
#include "plot/vtkplotClass.h"

#include "finiteElements/ElementList.h"
#include "finiteElements/GenericFiniteElement.h"
#include "pointercollection/pointercollection.h"

namespace HierAMuS {

PlotControl::PlotControl() {
  m_plot = std::make_shared<vtkPlotInterface>();
  m_init = false;
  m_initMesh = false;
}

void PlotControl::initialize(PointerCollection &pointers)
{
  if(!m_init){
    m_init = true;
    auto infos = pointers.getInfoData();
    m_plot->initFileNames(infos->fileNames[FileHandling::directory],
                          infos->fileNames[FileHandling::infile]);
  }
}

void PlotControl::initializeMesh(PointerCollection &pointers)
{
  if(!m_initMesh){
    auto elemList = pointers.getElementList();
    indexType numberOfElements = elemList->getNumberOfElements();

    for (auto i = 0; i < numberOfElements; ++i) {
      elemList->getElement(pointers, i)->toParaviewAdapter(pointers, *m_plot,
                                                 ParaviewSwitch::Mesh);
    }
  }
}



void PlotControl::timeUpdate(prec time)
{
  m_plot->TimeUpdate(time);
}

void PlotControl::toFile(PointerCollection &pointers) {

  auto infos = pointers.getInfoData();

  this->initialize(pointers);
  this->initializeMesh(pointers);

  auto elemList = pointers.getElementList();
  indexType numberOfElements = elemList->getNumberOfElements();

  for (auto i = 0; i < numberOfElements; ++i) {
    elemList->getElement(pointers, i)->toParaviewAdapter(pointers, *m_plot,
                                               ParaviewSwitch::Weights);
    elemList->getElement(pointers, i)
        ->toParaviewAdapter(pointers, *m_plot,
                                               ParaviewSwitch::Solution);
    elemList->getElement(pointers, i)
        ->toParaviewAdapter(pointers, *m_plot,
                                               ParaviewSwitch::ProjectedValues);
  }
  m_plot->toFile(infos->fileNames[FileHandling::directory],
                 infos->fileNames[FileHandling::infile]);
}
} // namespace HierAMuS