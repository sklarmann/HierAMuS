// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "control/ParameterList.h"
#include "datatypes.h"
#include "shapefunctions/LegendreShapes.h"
#include "solver/SolutionTypes.h"
#include "solver/StaticSolutionStateHomogenization.h"

#include <memory>
#include <pointercollection/pointercollection.h>

#include "EquationHandler.h"
#include <finiteElements/ElementList.h>
#include <geometry/GeometryData.h>
#include <materials/ElementformulationList.h>
#include <materials/MaterialList.h>
#include <materials/MaterialformulationList.h>
#include <plot/vtkplotClass.h>
#include "plot/plotControl.h"

#include <solver/GenericSolutionState.h>
#include <solver/StaticSolutionState.h>

#include <control/HandlingStructs.h>
#include "LoadList.h"
#include "PropfunctionHandler.h"

#include <algorithm>
#include <cstdlib>

#include <solver/StaticSolutionState.h>
#include <solver/TransientSolutionNewmark.h>

#include "control/BinaryWrite.h"

#include <filesystem>


#ifdef OPENMP
#include <omp.h>
#endif

namespace HierAMuS {

PointerCollection::PointerCollection() {
#ifdef OPENMP
  Eigen::initParallel();
#endif

  this->geometryData = nullptr;
  this->solutionState = nullptr;
  this->eqHandler = nullptr;
  this->elementList = nullptr;
  this->materialList = nullptr;
  this->props = nullptr;
  this->loads = nullptr;
  this->prescribedDisplacements = nullptr;
  this->vtkPlot = nullptr;
  this->plotControl = nullptr;
  this->materialFormulationList = nullptr;
  this->elementFormulationList = nullptr;

  InfoData temp;
  this->m_hasRVEs = false;

}

PointerCollection::~PointerCollection() {}

void PointerCollection::renew() {

  this->geometryData = nullptr;
  this->solutionState = nullptr;
  this->eqHandler = nullptr;
  this->elementList = nullptr;
  this->materialList = nullptr;
  this->props = nullptr;
  this->loads = nullptr;
  this->vtkPlot = nullptr;
  this->plotControl = nullptr;
  this->materialFormulationList = nullptr;
  this->elementFormulationList = nullptr;
  this->computationData = nullptr;
  this->newGeometry();
  this->newEqHandler();
  this->newElementList();
  this->newMaterialList();
  this->newLoadList();
  this->newPrescribedDisplacements();
  this->newVtkPlot();
  this->newMaterialFormulationList();
  this->newElementFormulationList();
  this->newComputationData();
}

void PointerCollection::setMaxThreads(indexType maxThreads) {
#ifdef OPENMP
  omp_set_num_threads(maxThreads);
#endif
}

void PointerCollection::newGeometry() {
  this->geometryData = std::make_shared<Geometry::GeometryData>();
  // if (this->geometryData == 0) {
  //  this->geometryData = new Geometry::GeometryData(this);
  //} else {
  //}
}

void PointerCollection::newEqHandler() {
  if (this->eqHandler == nullptr) {
    this->eqHandler = std::make_shared<EquationHandler>();
  } else {
  }
}

void PointerCollection::newElementList() {
  this->elementList = std::make_shared<FiniteElement::ElementList>();
}

void PointerCollection::newMaterialList() {
  this->materialList = std::make_shared<Materials::MaterialList>();
}


std::shared_ptr<PropfunctionHandler> PointerCollection::getPropLoads() {
  return this->solutionState->getProps();
}

void PointerCollection::newLoadList() {
  this->loads = std::make_shared<LoadList>();
}

void PointerCollection::newVtkPlot() {
  this->vtkPlot = std::make_shared<vtkPlotInterface>();
  this->plotControl = std::make_shared<PlotControl>();
}

void PointerCollection::newMaterialFormulationList() {
  this->materialFormulationList =
      std::make_shared<Materials::MaterialFormulationList>();
}

void PointerCollection::newElementFormulationList() {
  this->elementFormulationList =
      std::make_shared<Materials::ElementFormulationList>();
}

void PointerCollection::setInfoData(std::shared_ptr<InfoData> Info) {
  this->Info = Info;
}
std::shared_ptr<GenericSolutionState> PointerCollection::getSolutionState() {
  return this->solutionState;
}
void PointerCollection::solutionStateToFile(std::string fileName) {
  std::string dir =
      this->Info->fileNames[FileHandling::directory] + "solution/";

  if (!std::filesystem::exists(dir)) {
    std::filesystem::create_directory(dir);
  }

  std::string filedir = dir + fileName;
  std::ofstream out(filedir.c_str(), std::fstream::out | std::fstream::binary |
                                         std::fstream::trunc);

  if (!out.is_open())
    throw std::runtime_error("Could not open file " + filedir);

  this->getSPDLogger().debug("Writing solution to file " + filedir);

  auto soltype = this->solutionState->getType();
  out.write(reinterpret_cast<char *>(&soltype), sizeof(soltype));
  this->solutionState->toFile(*this, out);
  out.close();
}
void PointerCollection::solutionStateFromFile(std::string fileName) {

  std::string dir = this->Info->fileNames[FileHandling::directory];

  std::string filedir = dir + "solution/" + fileName;
  std::ifstream in(filedir.c_str(), std::ios::in | std::ios::binary);
  SolutionTypes soltype;
  this->getSPDLogger().debug("Reading solution from file " + filedir);
  if (!in.is_open())
    throw std::runtime_error("Could not open file " + filedir);
  in.read(reinterpret_cast<char *>(&soltype), sizeof(soltype));

  ParameterList dummyParam;
  switch (soltype) {
  case SolutionTypes::StaticSolutionState:
    this->getSPDLogger().debug("Reading static solution state");
    this->solutionState = std::make_shared<StaticSolutionState>(dummyParam);
    break;
  case SolutionTypes::StaticSolutionHomogenization:
    this->getSPDLogger().debug("Reading static homogenization solution state");
    this->solutionState =
        std::make_shared<StaticSolutionStateHomogenization>(dummyParam);

    break;
  default:
    throw std::runtime_error(
        "Uknown solution type when trying to read it from file");
  }
  this->solutionState->setProps(this->props);
  this->solutionState->fromFile(*this, in);

  in.close();
}
void PointerCollection::RVEDataToFile(std::string fileName)
{
  std::string dir =
      this->Info->fileNames[FileHandling::directory] + "solution/";

  if (!std::filesystem::exists(dir)) {
    std::filesystem::create_directory(dir);
  }

  std::string filedir = dir + fileName;
  std::ofstream out(filedir.c_str(), std::fstream::out | std::fstream::binary |
                                         std::fstream::trunc);

  if (!out.is_open())
    throw std::runtime_error("Could not open file " + filedir);

  this->getSPDLogger().debug("Writing solution to file " + filedir);
  this->solutionState->RVEDatatoFile(*this, out);
  out.close();
}
void PointerCollection::RVEDataFromFile(std::string fileName)
{
  std::string dir = this->Info->fileNames[FileHandling::directory];

  std::string filedir = dir + "solution/" + fileName;
  std::ifstream in(filedir.c_str(), std::ios::in | std::ios::binary);
  
  this->getSPDLogger().debug("Reading solution from file " + filedir);
  if (!in.is_open())
    throw std::runtime_error("Could not open file " + filedir);
 
  this->solutionState->RVEDatafromFile(*this, in);

  in.close();
}
void PointerCollection::setSolutionState(
    std::shared_ptr<GenericSolutionState> solutionStateIn) {
  if (!this->props) {
    this->props = std::make_shared<PropfunctionHandler>();
  }
  this->solutionState = solutionStateIn;
  this->solutionState->setProps(this->props);
}
std::shared_ptr<Geometry::GeometryData> PointerCollection::getGeometryData() {
  return this->geometryData;
}
void PointerCollection::setGeometryData(
    std::shared_ptr<Geometry::GeometryData> geoDataIn) {
  this->geometryData = geoDataIn;
}
void PointerCollection::setEquationHandler(
    std::shared_ptr<EquationHandler> eqHandler) {
  this->eqHandler = eqHandler;
}
std::shared_ptr<EquationHandler> PointerCollection::getEquationHandler() {
  return this->eqHandler;
}
void PointerCollection::setElementList(
    std::shared_ptr<FiniteElement::ElementList> elementList) {
  this->elementList = elementList;
}
std::shared_ptr<FiniteElement::ElementList>
PointerCollection::getElementList() {
  return this->elementList;
}
void PointerCollection::setMaterialList(
    std::shared_ptr<Materials::MaterialList> matList) {
  this->materialList = matList;
}
std::shared_ptr<Materials::MaterialList> PointerCollection::getMaterialList() {
  return this->materialList;
}
std::shared_ptr<InfoData> PointerCollection::getInfoData() {
  return this->Info;
}
std::shared_ptr<LoadList> PointerCollection::getLoadList() {
  return this->loads;
}

void PointerCollection::newPrescribedDisplacements()
{
  this->prescribedDisplacements = std::make_shared<LoadList>();
}

std::shared_ptr<LoadList> PointerCollection::getPrescribedDisplacements() {
  return this->prescribedDisplacements;
}

std::shared_ptr<vtkPlotInterface> PointerCollection::getVtkPlotInterface() {
  return this->vtkPlot;
}

auto PointerCollection::getPlotControlInterface()
    -> std::shared_ptr<PlotControl> {
  return this->plotControl;
}
void PointerCollection::newComputationData() {
  this->computationData = std::make_shared<ParameterList>();
}
void PointerCollection::setExternalComputationData(
    std::shared_ptr<ParameterList> compDataExt) {
  this->computationData = compDataExt;
}
auto PointerCollection::getCompuationData() -> std::shared_ptr<ParameterList> {
  return this->computationData;
}
void PointerCollection::setMaterialFormulationsList(
    std::shared_ptr<Materials::MaterialFormulationList> matlist) {
  this->materialFormulationList = matlist;
}
std::shared_ptr<Materials::MaterialFormulationList>
PointerCollection::getMaterialFormulationList() {
  return this->materialFormulationList;
}
void PointerCollection::setElementFormulationList(
    std::shared_ptr<Materials::ElementFormulationList> elemList) {
  this->elementFormulationList = elemList;
}
std::shared_ptr<Materials::ElementFormulationList>
PointerCollection::getElementFormulationList() {
  return this->elementFormulationList;
}



//auto PointerCollection::getLogger() -> OutputHandler & {
//  return this->Info->Log;
//}

auto PointerCollection::getSPDLogger() -> spdlog::logger & {
  return this->Info->Log.getSPDLogger();
}

auto PointerCollection::getShallowCopy() -> std::shared_ptr<PointerCollection> {
  std::shared_ptr<PointerCollection> newPointers = std::make_shared<PointerCollection>();
  newPointers->geometryData = this->geometryData;

  newPointers->eqHandler = this->eqHandler;
  newPointers->elementList = this->elementList;
  newPointers->loads = this->loads;
  newPointers->prescribedDisplacements = this->prescribedDisplacements;
  newPointers->Info = this->Info;
  newPointers->props = this->props;
  newPointers->plotControl = this->plotControl;
  newPointers->solutionState = this->solutionState->getCopy();

  newPointers->materialList = this->materialList;
  newPointers->materialFormulationList = this->materialFormulationList;
  newPointers->elementFormulationList = this->elementFormulationList;

  return newPointers;
}
void PointerCollection::setRVEs() {
  this->solutionState->setHasRVE();
  this->m_hasRVEs = true;
}
auto PointerCollection::hasRVEs() -> bool { return this->m_hasRVEs; };
void PointerCollection::setSolutionState(SolutionTypes type,
                                         ParameterList &paramList) {

  if (!this->props)
  {
    this->props = std::make_shared<PropfunctionHandler>();
  }
  switch (type) {
  case SolutionTypes::StaticSolutionState:
    this->solutionState = std::make_shared<StaticSolutionState>(paramList);
    this->solutionState->setProps(this->props);
    break;
  case SolutionTypes::StaticSolutionHomogenization:
    this->solutionState =
        std::make_shared<StaticSolutionStateHomogenization>(paramList);
    this->solutionState->setProps(this->props);
  case SolutionTypes::TransientSolutionNewmark:
    this->solutionState = std::make_shared<TransientSolutionNewmark>(paramList);
    this->solutionState->setProps(this->props);
    break;
  case HierAMuS::SolutionTypes::GenericSolutionState:
  case HierAMuS::SolutionTypes::LinearStaticSolutionState:
  case HierAMuS::SolutionTypes::Transient: {
    std::cout << "Selected Solution Type not implemented!" << std::endl;
  }
  }
}
void PointerCollection::setExternalVtkPlot(
    std::shared_ptr<vtkPlotInterface> vtkPlotExt) {
  this->vtkPlot = vtkPlotExt;
}




} /* namespace HierAMuS */
