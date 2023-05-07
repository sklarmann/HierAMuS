// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <vector>
#include <string>

#include <algorithm>
#include <sstream>

#include <plot/vtkplotClass.h>



#ifdef USE_VTK
#include <vtkActor.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
//#include <vtkExtractEdges.h>
#include <vtkFloatArray.h>
#include <vtkLine.h>
#include <vtkLookupTable.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
//#include <vtkUnstructuredGridGeometryFilter.h>
#include <vtkVersion.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#endif

#include <regex>

#include "catalystcxxinterface/management.h"

namespace HierAMuS {


vtkPlotInterface::vtkPlotInterface() {

  this->step = 0;

  this->adapterInit = false;
  this->paraViewManager = std::make_shared<managementClass>();
}




vtkUnstructuredGrid *vtkPlotInterface::getGrid() {
  return NULL;
}


inline void vtkPlotInterface::interact() {

}

void vtkPlotInterface::initFileNames(std::string folder, std::string file)
{
  std::string outputFile, pvdFile, tempFileName, adapterFileName;
  outputFile = folder;
  pvdFile = outputFile + "parvout/";
  outputFile += "parvout/model/";
  tempFileName = file;
  adapterFileName = folder + file;
  std::size_t pos = tempFileName.find(".txt");
  if (pos != std::string::npos) {
    tempFileName = tempFileName.substr(0, pos);
  }
  outputFile += tempFileName;
  pvdFile += tempFileName;
  std::stringstream tempStr;
  tempStr << this->step;
  ++this->step;
  outputFile += tempStr.str() + ".vtm";
  pvdFile += ".pvd";

  if(!this->adapterInit){
    this->adapterInit = true;
    char *c = new char[adapterFileName.size() + 1];
    strcpy(c, adapterFileName.c_str());
#ifdef USE_VTK
    this->paraViewManager->initialize(c);
#endif
    delete[] c;
  }
}


void vtkPlotInterface::toFile(
    std::string folder, std::string file) {

#ifdef USE_VTK
  vtkSmartPointer<vtkMultiBlockDataSet> block =
      vtkSmartPointer<vtkMultiBlockDataSet>::New();
#endif
  std::string outputFile, pvdFile, tempFileName, adapterFileName;
  outputFile = folder;
  pvdFile = outputFile + "parvout/";
  outputFile += "parvout/model/";
  tempFileName = file;
  adapterFileName = folder + file;
  std::size_t pos = tempFileName.find(".txt");
  if (pos != std::string::npos) {
    tempFileName = tempFileName.substr(0, pos);
  }
  outputFile += tempFileName;
  pvdFile += tempFileName;
  std::stringstream tempStr;
  tempStr << this->step;
  ++this->step;
  outputFile += tempStr.str() + ".vtm";
  pvdFile += ".pvd";

#ifdef USE_VTK
  this->paraViewManager->normalizeFieldByWeight();
  this->paraViewManager->writeFile();
#endif
}


void vtkPlotInterface::pvdFileReader(
    std::string &pvdFile, std::vector<std::string> &FNames,
    std::vector<prec> &timesteps) {
  FNames.clear();
  timesteps.clear();

  std::ifstream file(pvdFile);
  if (file.good()) {
    // file.open(fileName.str());

    std::vector<std::string> fcon;
    std::string line;

    while (std::getline(file, line)) {
      fcon.push_back(line);
    }

    int start = 3;
    int end = static_cast<int>(fcon.size()) - 3;

    std::regex time("timestep=\"([^\"]*)\"");
    std::regex fname("file=\"([^\"]*)\"");
    std::smatch st, sf;
    while (start < end) {
      line = fcon[start];
      if (std::regex_search(line, st, time)) {
        double time = std::stod(st[1]);
        line = fcon[start + 1];
        std::regex_search(line, sf, fname);
        timesteps.push_back(time);
        FNames.push_back(sf[1]);
      }

      start += 2;
    }
    if (this->step == 1) {
      FNames.erase(FNames.begin() + this->step - 1, FNames.end());
      timesteps.erase(timesteps.begin() + this->step - 1, timesteps.end());
    }
  }
  file.close();
}


void vtkPlotInterface::pvdFileWriter(
    std::string &pvdFile, std::vector<std::string> &FNames,
    std::vector<prec> &timesteps) {
  std::ofstream file;
  file.open(pvdFile);

  std::string space;
  space = "  ";

  file << "<?xml version =\"1.0\"?>\n"
       << "<VTKFile type=\"Collection\" version=\"0.1\">\n";

  file << space << "<Collection>\n";
  space = "    ";

  for (auto i = 0; i < this->step; ++i) {
    file << space << "<DataSet timestep=\"" << timesteps[i]
         << "\" group=\"\" part=\"\"\n"
         << space << "file=\"" << FNames[i]
         << "\"/>\n"; //.substr(pos+1,FNames[i].length())
  }

  space = "  ";
  file << space << "</Collection>\n</VTKFile>" << std::endl;

  file.close();
}


void vtkPlotInterface::timeUpdate(const prec &time) {
#ifdef USE_VTK
  double tt = static_cast<double>(time);
  this->paraViewManager->TimeUpdate(tt);
#endif
}






void vtkPlotInterface::initialize(std::string fileName) {
#ifdef USE_VTK
  auto temp = fileName.c_str();
  this->paraViewManager->initialize(temp);
#endif
}


void vtkPlotInterface::CoProcess() {
#ifdef USE_VTK
  this->paraViewManager->CoProcess();
#endif
}


void vtkPlotInterface::finalizeCoProcessing() {
#ifdef USE_VTK
  this->paraViewManager->finalizeCoProcessing();
#endif
}


void vtkPlotInterface::isRVEfunc() {
#ifdef USE_VTK
  this->paraViewManager->isRVEfunc();
#endif
}


void vtkPlotInterface::TimeUpdate(prec time) {
#ifdef USE_VTK
  double tt = static_cast<double>(time);
  this->paraViewManager->TimeUpdate(tt);
#endif
}


void vtkPlotInterface::writeFile() {
#ifdef USE_VTK
  this->paraViewManager->writeFile();
#endif
}


void vtkPlotInterface::addPoint(indexType main, indexType part, indexType id, prec x, prec y, prec z) {
#ifdef USE_VTK
  int cmain = static_cast<int>(main);
  int cpart = static_cast<int>(part);
  int cid = static_cast<int>(id);
  double cx = static_cast<double>(x);
  double cy = static_cast<double>(y);
  double cz = static_cast<double>(z);
  this->paraViewManager->addPoint(cmain, cpart, cid, cx, cy, cz);
#endif
}


void vtkPlotInterface::addPoint(indexType main, indexType part, indexType id,
                Types::Vector3<prec> &coors) {
#ifdef USE_VTK
  int cmain = static_cast<int>(main);
  int cpart = static_cast<int>(part);
  int cid = static_cast<int>(id);
  double cx = static_cast<double>(coors(0));
  double cy = static_cast<double>(coors(1));
  double cz = static_cast<double>(coors(2));
  this->paraViewManager->addPoint(cmain, cpart, cid, cx, cy, cz);
#endif
}



void vtkPlotInterface::addCell(indexType main,
                                                indexType part,
                                                indexType FeapCellNumber,
                                                indexType FeapSubCellNumber,
                                                std::vector<indexType> FeapPoints,
                                                indexType numpts,
                                                int vtkNumber) {
#ifdef USE_VTK
  int cmain = static_cast<int>(main);
  int cpart = static_cast<int>(part);
  int cFeapCellNumber = static_cast<int>(FeapCellNumber);
  int cFeapSubCellNumber = static_cast<int>(FeapSubCellNumber);
  int cnumpts = static_cast<int>(numpts);
  std::vector<int> cFeapPoints;
  cFeapPoints.reserve(FeapPoints.size());
  for (auto &i : FeapPoints) {
    cFeapPoints.push_back(static_cast<int>(i));
  }
  this->paraViewManager->addCell(cmain,
                                cpart,
                                cFeapCellNumber,
                                cFeapSubCellNumber,
                                &cFeapPoints[0],
                                &cnumpts,
                                &vtkNumber);
#endif

}


inline void vtkPlotInterface::addCompleteFeapPointArray(indexType main,
                                                                         indexType part,
                                                                         std::vector<indexType> data,
                                                                         indexType num_comp,
                                                                         indexType numPoints,
                                                                         std::string name) {
#ifdef USE_VTK
  int cmain = static_cast<int>(main);
  int cpart = static_cast<int>(part);
  int cnum_comp = static_cast<int>(num_comp);
  int cnumPoints = static_cast<int>(numPoints);
  std::vector<int> cdata;
  cdata.reserve(data.size());
  for (auto &i : data) {
    cdata.push_back(static_cast<int>(i));
  }
  this->paraViewManager->addCompleteFeapPointArray(cmain, cpart, &cdata[0], cnum_comp, cnumPoints, name.c_str());
#endif
}


inline void vtkPlotInterface::addCompleteFeapPointArray(indexType main, indexType part,
                                                                         std::vector<prec> data, indexType num_comp,
                                                                         indexType numPoints, std::string name) {
#ifdef USE_VTK
  int cmain = static_cast<int>(main);
  int cpart = static_cast<int>(part);
  int cnum_comp = static_cast<int>(num_comp);
  int cnumPoints = static_cast<int>(numPoints);
  std::vector<double> cdata;
  cdata.reserve(data.size());
  for (auto &i : data) {
    cdata.push_back(static_cast<double>(i));
  }
  this->paraViewManager->addCompleteFeapPointArray(cmain, cpart, &cdata[0], cnum_comp, cnumPoints, name.c_str());
#endif
}


template<typename T>
void vtkPlotInterface::addCompleteFeapPointArrayToAll(indexType main,
                                                                       const std::vector<T> &data,
                                                                       indexType num_comp,
                                                                       indexType numPoints,
                                                                       const std::string &name) {
#ifdef USE_VTK
  int cmain = static_cast<int>(main);
  int cnum_comp = static_cast<int>(num_comp);
  int cnumPoints = static_cast<int>(numPoints);
  if (std::is_same<T, indexType>::value) {
    std::vector<int> cdata;
    cdata.reserve(data.size());
    for (auto &i : data) {
      cdata.push_back(static_cast<int>(i));
    }
    this->paraViewManager->addCompleteFeapPointArrayToAll(cmain, &cdata[0], cnum_comp, cnumPoints, name.c_str());
  } else if (std::is_same<T, prec>::value) {
    std::vector<double> cdata;
    cdata.reserve(data.size());
    for (auto &i : data) {
      cdata.push_back(static_cast<double>(i));
    }
    this->paraViewManager->addCompleteFeapPointArrayToAll(cmain, &cdata[0], cnum_comp, cnumPoints, name.c_str());
  }

#endif
}


void vtkPlotInterface::setPointData(indexType main,
                                                     indexType part,
                                                     indexType feapId,
                                                     std::vector<indexType> data,
                                                     indexType num_comp,
                                                     std::string name) {
#ifdef USE_VTK
  int cmain = static_cast<int>(main);
  int cpart = static_cast<int>(part);
  int cpointId = static_cast<int>(feapId);
  int cnum_comp = static_cast<int>(num_comp);
  std::vector<int> cdata;
  cdata.reserve(data.size());
  for (auto &i : data) {
    cdata.push_back(static_cast<int>(i));
  }
  this->paraViewManager->setPointData(cmain, cpart, cpointId, &cdata[0], cnum_comp, name.c_str());
#endif
}


void vtkPlotInterface::setPointData(indexType main,
                                                     indexType part,
                                                     indexType feapId,
                                                     std::vector<prec> data,
                                                     indexType num_comp,
                                                     std::string name) {
#ifdef USE_VTK
  int cmain = static_cast<int>(main);
  int cpart = static_cast<int>(part);
  int cpointId = static_cast<int>(feapId);
  int cnum_comp = static_cast<int>(num_comp);
  std::vector<double> cdata;
  cdata.reserve(data.size());
  for (auto &i : data) {
    cdata.push_back(static_cast<double>(i));
  }
  this->paraViewManager->setPointData(cmain, cpart, cpointId, &cdata[0], cnum_comp, name.c_str());
#endif
}


void vtkPlotInterface::SumPointDataWeighted(indexType main,
                                                     indexType part,
                                                     std::vector<prec> &data,
                                                     indexType feapId,
                                                     indexType num_comp,
                                                     const std::string &name) {
#ifdef USE_VTK
  int cmain = static_cast<int>(main);
  int cpart = static_cast<int>(part);
  int cpointId = static_cast<int>(feapId);
  int cnum_comp = static_cast<int>(num_comp);
  std::vector<double> cdata;
  cdata.reserve(data.size());
  for (auto &i : data) {
    cdata.push_back(static_cast<double>(i));
  }
  this->paraViewManager->SumPointDataWeighted(cmain, cpart, &cdata[0], cpointId, cnum_comp, name.c_str());

#endif

}


void vtkPlotInterface::SumPointDataWeighted(indexType main,
                                                     indexType part,
                                                     std::vector<indexType> &data,
                                                     indexType feapId,
                                                     indexType num_comp,
                                                     const std::string &name) {
#ifdef USE_VTK
  int cmain = static_cast<int>(main);
  int cpart = static_cast<int>(part);
  int cpointId = static_cast<int>(feapId);
  int cnum_comp = static_cast<int>(num_comp);
  std::vector<int> cdata;
  cdata.reserve(data.size());
  for (auto &i : data) {
    cdata.push_back(static_cast<int>(i));
  }
  this->paraViewManager->SumPointDataWeighted(cmain, cpart, &cdata[0], cpointId, cnum_comp, name.c_str());

#endif
}



void vtkPlotInterface::SumPointData(indexType main,
                                                     indexType part,
                                                     std::vector<prec> &data,
                                                     indexType feapId,
                                                     indexType num_comp,
                                                     const std::string &name) {
#ifdef USE_VTK
  int cmain = static_cast<int>(main);
  int cpart = static_cast<int>(part);
  int cpointId = static_cast<int>(feapId);
  int cnum_comp = static_cast<int>(num_comp);
  std::vector<double> cdata;
  cdata.reserve(data.size());
  for (auto &i : data) {
    cdata.push_back(static_cast<double>(i));
  }
  this->paraViewManager->SumPointData(cmain, cpart, &cdata[0], cpointId, cnum_comp, name.c_str());

#endif

}


void vtkPlotInterface::SumPointData(indexType main,
                                                     indexType part,
                                                     std::vector<indexType> &data,
                                                     indexType feapId,
                                                     indexType num_comp,
                                                     const std::string &name) {
#ifdef USE_VTK
  int cmain = static_cast<int>(main);
  int cpart = static_cast<int>(part);
  int cpointId = static_cast<int>(feapId);
  int cnum_comp = static_cast<int>(num_comp);
  std::vector<int> cdata;
  cdata.reserve(data.size());
  for (auto &i : data) {
    cdata.push_back(static_cast<int>(i));
  }
  this->paraViewManager->SumPointData(cmain, cpart, &cdata[0], cpointId, cnum_comp, name.c_str());

#endif

}

} // namespace HierAMuS


