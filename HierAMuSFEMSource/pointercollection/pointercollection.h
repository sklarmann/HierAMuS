// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once
#include <Base/FEMBase.h>
#include <datatypes.h>
#include <forwarddeclaration.h>
#include <iostream>
#include "control/OutputHandler.h"
#include "plot/plotControl.h"
#include "plot/vtkplotClass.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "shapefunctions/IntegrationsPoints/helperClasses/IntegrationPointsManagement.h"
#include <shapefunctions/LegendreShapes.h>
#include <shapefunctions/LobattoShapes.h>
#include <solver/GenericSolutionState.h>
#include <solver/SolutionTypes.h>

#include <control/ParameterList.h>

#include <memory>
#include <spdlog/logger.h>


namespace HierAMuS {

class PlotControl;
class PointerCollection : public FEMBase {
public:
  PointerCollection();
  virtual ~PointerCollection();

  virtual void renew();

  virtual void setMaxThreads(indexType maxThreads);

  virtual void
  setGeometryData(std::shared_ptr<Geometry::GeometryData> geoDataIn);
  virtual std::shared_ptr<Geometry::GeometryData> getGeometryData();
  virtual void newGeometry();

  virtual void
  setSolutionState(std::shared_ptr<GenericSolutionState> solutionStateIn);
  virtual void setSolutionState(SolutionTypes type, ParameterList &paramList);
  virtual std::shared_ptr<GenericSolutionState> getSolutionState();
  virtual void solutionStateToFile(std::string fileName);
  virtual void solutionStateFromFile(std::string fileName);

  
  virtual void RVEDataToFile(std::string fileName);
  virtual void RVEDataFromFile(std::string fileName);
  // void newSolutionState();

  virtual void setEquationHandler(std::shared_ptr<EquationHandler> eqHandler);
  virtual std::shared_ptr<EquationHandler> getEquationHandler();
  virtual void newEqHandler();

  virtual void
  setElementList(std::shared_ptr<FiniteElement::ElementList> elementList);
  virtual std::shared_ptr<FiniteElement::ElementList> getElementList();
  virtual void newElementList();

  virtual void
  setMaterialList(std::shared_ptr<Materials::MaterialList> matList);
  virtual std::shared_ptr<Materials::MaterialList> getMaterialList();
  virtual void newMaterialList();


  virtual void setInfoData(std::shared_ptr<InfoData> Info);
  virtual std::shared_ptr<InfoData> getInfoData();

  virtual std::shared_ptr<PropfunctionHandler> getPropLoads();

  virtual void newLoadList();
  virtual std::shared_ptr<loadList> getLoadList();

  virtual void newPrescribedDisplacements();
  virtual std::shared_ptr<loadList> getPrescribedDisplacements();

  virtual void newVtkPlot();
  virtual void setExternalVtkPlot(std::shared_ptr<vtkPlotInterface> vtkPlotExt);
  virtual auto getVtkPlotInterface() -> std::shared_ptr<vtkPlotInterface>;
  virtual auto getPlotControlInterface() -> std::shared_ptr<PlotControl>;

  virtual void newComputationData();
  virtual void
  setExternalComputationData(std::shared_ptr<ParameterList> compDataExt);
  virtual auto getCompuationData() -> std::shared_ptr<ParameterList>;

  virtual void setMaterialFormulationsList(
      std::shared_ptr<Materials::MaterialFormulationList> matlist);
  virtual void newMaterialFormulationList();
  virtual auto getMaterialFormulationList()
      -> std::shared_ptr<Materials::MaterialFormulationList>;

  virtual void setElementFormulationList(
      std::shared_ptr<Materials::ElementFormulationList> elemList);
  virtual void newElementFormulationList();
  virtual auto getElementFormulationList()
      -> std::shared_ptr<Materials::ElementFormulationList>;

  //static IntegrationPoints getIntegrationPoints();
  static IntegrationPoints getIntegrationPoints(indexType elementId);
  static LobattoShapes &getLobattoShapes();
  static LegendreShapes &getLegendreShapes();

  // auto getLogger() -> OutputHandler &;

  auto getSPDLogger() -> spdlog::logger &;

  auto getShallowCopy() -> std::shared_ptr<PointerCollection>;

  void setRVEs()
  {
    this->solutionState->setHasRVE();
    this->m_hasRVEs = true;
  };
  auto hasRVEs() -> bool { return this->m_hasRVEs; };

private:
  std::shared_ptr<Geometry::GeometryData> geometryData;
  std::shared_ptr<GenericSolutionState> solutionState;
  std::shared_ptr<EquationHandler> eqHandler;
  std::shared_ptr<FiniteElement::ElementList> elementList;
  std::shared_ptr<loadList> loads;
  std::shared_ptr<loadList> prescribedDisplacements;
  std::shared_ptr<InfoData> Info;
  std::shared_ptr<PropfunctionHandler> props;
  std::shared_ptr<vtkPlotInterface> vtkPlot;
  std::shared_ptr<PlotControl> plotControl;
  std::shared_ptr<ParameterList> computationData;

  std::shared_ptr<Materials::MaterialList> materialList;
  std::shared_ptr<Materials::MaterialFormulationList> materialFormulationList;
  std::shared_ptr<Materials::ElementFormulationList> elementFormulationList;


  bool m_hasRVEs;

  static IntegrationPointsManagement intPoints;
  static LobattoShapes LobattoFunctions;
  static LegendreShapes LegendreFunctions;

};

} /* namespace HierAMuS */
