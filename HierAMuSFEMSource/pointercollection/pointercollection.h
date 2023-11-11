// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once
#include "datatypes.h"

#include "shapefunctions/IntegrationsPoints/helperClasses/IntegrationPointsManagement.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "solver/SolutionTypes.h"

#include <memory>
#include <iostream>



namespace spdlog {
  class logger;
}

namespace HierAMuS {

class EquationHandler;
class LoadList;
class PlotControl;
class ParameterList;
class PropfunctionHandler;
class vtkPlotInterface;
class GenericSolutionState;

namespace Geometry {
class GeometryData;
}

namespace FiniteElement {
class ElementList;
}

namespace Materials {
class ElementFormulationList;
class MaterialList;
class MaterialFormulationList;
}

class PointerCollection {
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
  virtual std::shared_ptr<LoadList> getLoadList();

  virtual void newPrescribedDisplacements();
  virtual std::shared_ptr<LoadList> getPrescribedDisplacements();

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


  // auto getLogger() -> OutputHandler &;

  auto getSPDLogger() -> spdlog::logger &;

  auto getShallowCopy() -> std::shared_ptr<PointerCollection>;

  void setRVEs();
  auto hasRVEs() -> bool;

private:
  std::shared_ptr<Geometry::GeometryData> geometryData;
  std::shared_ptr<GenericSolutionState> solutionState;
  std::shared_ptr<EquationHandler> eqHandler;
  std::shared_ptr<FiniteElement::ElementList> elementList;
  std::shared_ptr<LoadList> loads;
  std::shared_ptr<LoadList> prescribedDisplacements;
  std::shared_ptr<InfoData> Info;
  std::shared_ptr<PropfunctionHandler> props;
  std::shared_ptr<vtkPlotInterface> vtkPlot;
  std::shared_ptr<PlotControl> plotControl;
  std::shared_ptr<ParameterList> computationData;

  std::shared_ptr<Materials::MaterialList> materialList;
  std::shared_ptr<Materials::MaterialFormulationList> materialFormulationList;
  std::shared_ptr<Materials::ElementFormulationList> elementFormulationList;


  bool m_hasRVEs;


};

} /* namespace HierAMuS */
