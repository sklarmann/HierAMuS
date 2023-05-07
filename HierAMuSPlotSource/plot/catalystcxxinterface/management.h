// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>
#include <vector>

#ifdef USE_VTK
#include <vector>
#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>


#include <plot/catalystcxxinterface/gridHandler.h>
#include <vtkSMProxy.h>
#include <vtkSMSourceProxy.h>


#include <vtkXMLMultiBlockDataWriter.h>

#include <regex>
#include <vtkXMLPVAnimationWriter.h>
#include <vtkXMLPVDWriter.h>

#include <plot/catalystcxxinterface/gridHandler.h>
#include <vtkCPProcessor.h>
#include <vtkInformation.h>
#include <vtkLiveInsituLink.h>
#include <vtkPVTrivialProducer.h>
#include <vtkSMInputProperty.h>
#include <vtkSMProxyManager.h>
#include <vtkSMSessionProxyManager.h>

#include <vtkTable.h>
#endif
class managementClass {
public:
  managementClass();
  ~managementClass();

  // Deprecated

  void addPoints(double *coor, int &numpoints, int &ndm);

  vtkUnstructuredGrid *getGrid(const int &main, const int &part);

  //

  void initialize(const char *sname);
  void CoProcess();
  void finalizeCoProcessing();
  void isRVEfunc();

  void TimeUpdate(double &time);

  void writeFile();

  void addPoint(const int &main, const int &part, const int &id,
                const double &x, const double &y, const double &z);
  void addCell(const int &main, const int &part, int &FeapCellNumber,
               int &FeapSubCellNumber, int *FeapPoints, int *numpts,
               int *vtkNumber);

  template <typename T>
  void addCompleteFeapPointArray(const int &main, const int &part,
                                 const T *data, const int &num_comp,
                                 const int &numPoints, const char *name);
  template <typename T>
  void addCompleteFeapPointArrayToAll(const int &main, const T *data,
                                      const int &num_comp, const int &numPoints,
                                      const char *name);
  template <typename T>
  void setPointData(const int &main, const int &part, const int &feapId,
                    const T *data, const int &num_comp, const char *name);
  template <typename T>
  void SumPointDataWeighted(const int &main, const int &part, const T *data,
                    const int &feapId, const int &num_comp, const char *name);
  template <typename T>
  void SumPointData(const int &main, const int &part, const T *data,
                    const int &feapId, const int &num_comp, const char *name);
  template <typename T>
  void setCellData(const int &main, const int &part, const int &feapCellNumber,
                   const int &feapSubCellNumber, const T *data,
                   const int &num_comp, const char *name);

  template <typename T>
  void SetFieldData(const int &main, const T *data, const int &num_comp,
                    const char *name);

  void Update();

  void SetAllToZero();

  void normalizeFieldByWeight();

private:

  bool checkGrid(const int &main, const int &part);
  void addGrid(const int &main, const int &part);
  void checkTploArray(const char *name);
  void checkTploSize();
  void addTploValue(const char *name, double &value);
#ifdef USE_VTK
  vtkSmartPointer<vtkMultiBlockDataSet> m_toplot;
  std::map<int, std::map<int, gridHandler>> gridMapper;
#endif
  // int m_numGrids;

  // bool m_connected;

  bool isRVE;
#ifdef USE_VTK
  bool coProcInit;
  // bool m_meshPresent;
  // bool m_doCoProcessing;
  // bool m_dataChanged;

  double m_time;
  int m_tstep, m_needCoProc;

  bool m_paraIsInit;
  vtkCPProcessor *m_proc;
  vtkLiveInsituLink *m_link;
  vtkSMSessionProxyManager *m_spxm;
  vtkSMProxy *m_px;
  vtkSMSourceProxy *m_spx;

  vtkSMProxy *m_px2;
  vtkSMSourceProxy *m_spx2;

  vtkPVTrivialProducer *m_meshProd;
  vtkPVTrivialProducer *m_tploProd;

  vtkSmartPointer<vtkTable> m_tplo;
  vtkIdType m_row;
#endif
  /// File Writer
  void pvdCollectionWriter();
  void pvdCollectionReader();
  std::string m_baseFileName, m_path;
  std::map<double, std::string> writerTimeMap;
  bool fileWritten;
  int writerStep;

};

/**
        Definition of Contructor and Destructor
                Constructor:
                                                - Initializes the
vtkMultiBlockDataSet which contains the different meshes
                                                - Initializes the following
variables: coProcInit = false:		Initially coprocessing is inactive
                                                                        m_tstep
= 0:			Current timestep m_time = 0.0: Current time isRVE =
false:			Wheter plot-data represent a RVE

**/
inline managementClass::managementClass() {
#ifdef USE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> grid =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

  this->m_toplot = vtkSmartPointer<vtkMultiBlockDataSet>::New();

  this->coProcInit = false;
  // this->m_meshPresent = false;

  this->m_tstep = 0;
  this->m_time = 0.0;
  this->isRVE = false;

  this->m_tplo = vtkSmartPointer<vtkTable>::New();
  this->m_row = 0;
  this->m_proc = 0;
  this->writerStep = 0;
  this->m_needCoProc = 0;
  this->m_paraIsInit = false;

  this->m_spx = 0;
  this->m_spx2 = 0;
  this->m_px = 0;
  this->m_px2 = 0;
  this->m_spxm = 0;
  this->m_link = 0;
  this->fileWritten = false;

  this->m_meshProd = 0;
  this->m_tploProd = 0;
#endif
}

/**
                Destructor: Currently nothing to do here.
**/
inline managementClass::~managementClass() {
  // this->writer->Delete();
}

/**
                initialize clears the data, reads and sets up the file name.
**/
inline void managementClass::initialize(const char *sname) {
#ifdef USE_VTK
  // File writing par
  std::string temp;
  temp = sname;

  std::replace(temp.begin(), temp.end(), '\\', '/');
  // std::cout << temp;

  auto pos = temp.find_last_of('/');
  if (pos != std::string::npos) {
    this->m_path = temp.substr(0, pos) + "/";
    this->m_baseFileName = temp.substr(pos + 1, temp.length());
  } else {
    this->m_path = "";
    this->m_baseFileName = temp;
  }

  // pos = this->m_baseFileName.find_last_of('.');
  // this->m_baseFileName = this->m_baseFileName.substr(0, pos);

  this->m_baseFileName[0] = 'R';
  // std::cout << "\n" << this->m_path << "\n" << this->m_baseFileName <<
  // std::endl;

  this->fileWritten = false;
  this->writerStep = 0;
  this->writerTimeMap.clear();
  this->m_time = 0.0;

  this->m_toplot->ReleaseData();
  for (auto it : this->gridMapper) {
    it.second.clear();
  }
  this->gridMapper.clear();

  this->m_tplo->ReleaseData();
  this->m_row = 0;

  char *buffer = new char[5];
  buffer[0] = 'T';
  buffer[1] = 'i';
  buffer[2] = 'm';
  buffer[3] = 'e';

  buffer[4] = '\0';
  ++this->m_row;
  this->addTploValue(buffer, this->m_time);
  // this->checkTploArray(buffer);
  // this->m_tplo->InsertNextBlankRow();
  //++this->m_row;
  // this->addTploValue(buffer, this->m_time);
  delete[] buffer;
#endif
}

inline void managementClass::SetAllToZero() {
#ifdef USE_VTK
  for (auto i = this->gridMapper.begin(); i != this->gridMapper.end(); ++i) {
    for (auto j = i->second.begin(); j != i->second.end(); ++j) {
      j->second.setAllToZero();
    }
  }
#endif
}

inline void managementClass::normalizeFieldByWeight() 
{
#ifdef USE_VTK
  for (auto i = this->gridMapper.begin(); i != this->gridMapper.end(); ++i) {
    for (auto j = i->second.begin(); j != i->second.end(); ++j) {
      j->second.normalizeFieldByWeight();
    }
  }
#endif
}

inline bool managementClass::checkGrid(const int &main, const int &part) {
#ifdef USE_VTK
  auto search = this->gridMapper.find(main);
  if (search != this->gridMapper.end()) {
    auto search2 = search->second.find(part);
    if (search2 != search->second.end())
      return true;
  }
#endif
  return false;
}

template <typename T>
inline void managementClass::addCompleteFeapPointArray(
    const int &main, const int &part, const T *data, const int &num_comp,
    const int &numPoints, const char *name) {
#ifdef USE_VTK
  if (this->checkGrid(main, part)) {
    this->gridMapper[main][part].CompleteFeapPointArray<T>(data, num_comp,
                                                           numPoints, name);
  } else {
    std::cout << "Error when trying to add complete data field " << name
              << " to main " << main << ", part " << part << "."
              << " Plot mesh does not exist!" << std::endl;
  }
#endif
  // this->Update();
}

template <typename T>
inline void managementClass::addCompleteFeapPointArrayToAll(
    const int &main, const T *data, const int &num_comp, const int &numPoints,
    const char *name) {
#ifdef USE_VTK
  auto it = this->gridMapper.find(main);
  if (it != this->gridMapper.end()) {
    for (auto &it2 : it->second) {
      it2.second.CompleteFeapPointArray<T>(data, num_comp, numPoints, name);
    }
  }
  this->Update();
#endif
}

template <typename T>
inline void managementClass::setPointData(const int &main, const int &part,
                                          const int &feapId, const T *data,
                                          const int &num_comp,
                                          const char *name) {
#ifdef USE_VTK
  /**
  int tlen = name_len + 1;
  char *buff = new char[name_len + 1];
  for (auto i = 0; i < name_len; ++i) {
          buff[i] = name[i];
  }
  buff[name_len] = '\0';
**/
  if (this->checkGrid(main, part)) {
    this->gridMapper[main][part].SetPointData<T>(data, feapId, num_comp, name);
  } else {
    std::cout << "Error when trying to add complete data field " << name
              << " to main " << main << ", part " << part << "."
              << " Plot mesh does not exist!" << std::endl;
  }

  this->Update();
#endif 
}

template <typename T>
inline void managementClass::SumPointDataWeighted(const int &main, const int &part,
                                          const T *data, const int &feapId,
                                          const int &num_comp,
                                          const char *name) {
#ifdef USE_VTK
  if (this->checkGrid(main, part)) {
    this->gridMapper[main][part].SumPointDataWeighted<T>(data, feapId, num_comp, name);
  } else {
    std::cout << "Error when trying to add complete data field " << name
              << " to main " << main << ", part " << part << "."
              << " Plot mesh does not exist!" << std::endl;
  }

  this->Update();
#endif
}

template <typename T>
inline void managementClass::SumPointData(const int &main, const int &part,
                                          const T *data, const int &feapId,
                                          const int &num_comp,
                                          const char *name) {
#ifdef USE_VTK
  if (this->checkGrid(main, part)) {
    this->gridMapper[main][part].SumPointData<T>(data, feapId, num_comp, name);
  } else {
    std::cout << "Error when trying to add complete data field " << name
              << " to main " << main << ", part " << part << "."
              << " Plot mesh does not exist!" << std::endl;
  }

  this->Update();
#endif
}

template <typename T>
inline void managementClass::SetFieldData(const int &main, const T *data,
                                          const int &num_comp,
                                          const char *name) {
#ifdef USE_VTK
  auto it = this->gridMapper.find(main);
  if (it != this->gridMapper.end()) {
    for (auto &it2 : it->second) {
      it2.second.SetFieldData(data, num_comp, name);
    }
  }
  // this->checkTploArray(name);
  double addData = static_cast<double>(data[0]);
  this->addTploValue(name, addData);
#endif
}

inline void managementClass::isRVEfunc() { this->isRVE = true; }

inline void managementClass::TimeUpdate(double &time) {
#ifdef USE_VTK
  ++this->m_tstep;
  this->m_time = time;
  if (!this->isRVE) {
    if (this->fileWritten) {
      ++this->writerStep;
      this->fileWritten = false;
    }
  }

  char *buffer = new char[5];
  buffer[0] = 'T';
  buffer[1] = 'i';
  buffer[2] = 'm';
  buffer[3] = 'e';

  buffer[4] = '\0';
  ++this->m_row;
  this->addTploValue(buffer, this->m_time);
  // this->checkTploArray(buffer);
  // this->m_tplo->InsertNextBlankRow();
  //++this->m_row;
  // this->addTploValue(buffer, this->m_time);
  delete[] buffer;
  this->SetAllToZero();
#endif
}

inline vtkUnstructuredGrid *managementClass::getGrid(const int &main,
                                                     const int &part) {
#ifdef USE_VTK
  this->checkGrid(main, part);
  vtkMultiBlockDataSet *mainBlock =
      vtkMultiBlockDataSet::SafeDownCast(this->m_toplot->GetBlock(main));
  vtkUnstructuredGrid *grid =
      vtkUnstructuredGrid::SafeDownCast(mainBlock->GetBlock(part));
  return grid;
#else
  return nullptr;
#endif
}

inline void managementClass::addPoint(const int &main, const int &part,
                                      const int &id, const double &x,
                                      const double &y, const double &z) {
#ifdef USE_VTK
  if (!this->checkGrid(main, part))
    this->addGrid(main, part);
  this->gridMapper[main][part].addPoint(id, x, y, z);
#endif
}

inline void managementClass::addCell(const int &main, const int &part,
                                     int &FeapCellNumber,
                                     int &FeapSubCellNumber, int *FeapPoints,
                                     int *numpts, int *vtkNumber) {
#ifdef USE_VTK
  if (this->checkGrid(main, part)) {
    this->gridMapper[main][part].addCell(FeapCellNumber, FeapSubCellNumber,
                                         FeapPoints, *numpts, *vtkNumber);
  } else {
    std::cout << "Error when trying to add cell to main " << main << ", part "
              << part << ". First initialize the mesh by adding its points."
              << std::endl;
  }
#endif
}

inline void managementClass::addPoints(double *coor, int &numpoints, int &ndm) {
#ifdef USE_VTK
  vtkUnstructuredGrid *grid =
      vtkUnstructuredGrid::SafeDownCast(this->m_toplot->GetBlock(0));

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  int k = 0;
  for (auto i = 0; i < numpoints; ++i) {
    std::cout << k << " " << coor[k] << " " << coor[k + 1] << " " << coor[k + 2]
              << std::endl;
    if (ndm == 3) {
      points->InsertNextPoint(coor[k], coor[k + 1], coor[k + 2]);
    } else if (ndm == 2) {
      points->InsertNextPoint(coor[k], coor[k + 1], 0.0);
    } else if (ndm == 1) {
      points->InsertNextPoint(coor[k], 0.0, 0.0);
    }
    k = k + ndm;
  }

  grid->SetPoints(points);
#endif
}

inline void managementClass::addGrid(const int &main, const int &part) {
#ifdef USE_VTK
  std::string blockname;

  if (main == 0) {
    std::stringstream temp;
    temp << "Material " << part;

    blockname = temp.str();
  }

  if (!this->checkGrid(main, part)) {
    auto search = this->gridMapper.find(main);
    if (search != this->gridMapper.end()) {

      vtkMultiBlockDataSet *block =
          vtkMultiBlockDataSet::SafeDownCast(this->m_toplot->GetBlock(main));
      block->SetBlock(part, this->gridMapper[main][part].getGrid());
      block->GetMetaData(part)->Set(vtkCompositeDataSet::NAME(),
                                    blockname.c_str());

    } else {
      this->gridMapper[main][part];

      this->m_toplot->SetBlock(main, vtkMultiBlockDataSet::New());
      vtkMultiBlockDataSet *block =
          vtkMultiBlockDataSet::SafeDownCast(this->m_toplot->GetBlock(main));
      block->SetBlock(part, this->gridMapper[main][part].getGrid());
      block->GetMetaData(part)->Set(vtkCompositeDataSet::NAME(),
                                    blockname.c_str());

      if (main == 0) {
        std::stringstream temp;
        temp << "Hauptmodell";
        this->m_toplot->GetMetaData(main)->Set(vtkCompositeDataSet::NAME(),
                                               temp.str().c_str());
      }
    }
  }
#endif
}

template <typename T>
inline void managementClass::setCellData(const int &main, const int &part,
                                         const int &feapCellNumber,
                                         const int &feapSubCellNumber,
                                         const T *data, const int &num_comp,
                                         const char *name) {
#ifdef USE_VTK
  if (this->checkGrid(main, part)) {
    this->gridMapper[main][part].SetCellData<T>(data, num_comp, feapCellNumber,
                                                feapSubCellNumber, name);
  } else {
    std::cout << "Error when trying to add complete data field " << name
              << " to main " << main << ", part " << part << "."
              << " Plot mesh does not exist!" << std::endl;
  }
#endif
}

/**
        File writing mehtods to store data on hard disk
**/
inline void managementClass::writeFile() {
#ifdef USE_VTK
  std::string wfilename;
  std::string folder;
  if (this->isRVE)
    this->pvdCollectionReader();

  std::stringstream temp, temp2, relName;
  relName << this->m_baseFileName << "/" << this->m_baseFileName
          << this->writerStep << ".vtm";
  temp << this->m_path << "parvout/" << this->m_baseFileName << "/";
  folder = temp.str();
  temp << this->m_baseFileName << this->writerStep << ".vtm";
  wfilename = temp.str();

  vtkSmartPointer<vtkXMLMultiBlockDataWriter> pwriter =
      vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
  pwriter->SetFileName(wfilename.c_str());
  pwriter->SetInputData(this->m_toplot);

  pwriter->SetEncodeAppendedData(1);
  pwriter->SetDataModeToAppended();

  // pwriter->SetDataModeToAscii();
  // pwriter->SetCompressorTypeToLZ4();
  pwriter->Write();

  this->writerTimeMap[this->m_time] = relName.str();
  this->fileWritten = true;

  this->pvdCollectionWriter();
#endif
}

/**
                pvdCollectionReader reads the pvd collection file and sets up
the corresponding data fields. This is only necessary in case of an RVE, because
data is lost due to restart of FEAP.
**/
inline void managementClass::pvdCollectionReader() {
#ifdef USE_VTK
  std::stringstream fileName;
  fileName << this->m_path << "parvout/";
  fileName << this->m_baseFileName << ".pvd";

  std::ifstream file(fileName.str());
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
        this->writerTimeMap[time] = sf[1];
        ++this->writerStep;
      }

      // line = fcon[start + 1];
      // std::regex_match(line, sf, fname);
      // for (auto i = 0; i < sf.size(); ++i) {
      //	std::cout << i << ": " << sf[i];
      //}
      start += 2;
    }
  }
  file.close();
#endif
}

/**
                pvdCollectionWriter writes the collection file which can be
loaded by Paraview.
**/
inline void managementClass::pvdCollectionWriter() {
#ifdef USE_VTK
  std::stringstream fileName;
  fileName << this->m_path << "parvout/";
  fileName << this->m_baseFileName << ".pvd";

  std::ofstream file;
  file.open(fileName.str());

  std::string space;
  space = "  ";

  file << "<?xml version =\"1.0\"?>\n"
       << "<VTKFile type=\"Collection\" version=\"0.1\">\n";

  file << space << "<Collection>\n";
  space = "    ";

  for (auto it : this->writerTimeMap) {
    file << space << "<DataSet timestep=\"" << it.first
         << "\" group=\"\" part=\"\"\n"
         << space << "file=\"" << it.second << "\"/>\n";
  }

  space = "  ";
  file << space << "</Collection>\n</VTKFile>" << std::endl;

  file.close();
#endif
}

/**
        CoProcessing Methods:
                CoProcess:
                                        - If link is not initialized, then this
routine generates the necessary objects to connect to Paraview. Current
implementation connects to "localhost" with port 22222.
                                        - After link is initialized, method
Update is called
**/
inline void managementClass::CoProcess() {
#ifdef USE_VTK
  if (!this->coProcInit) {
    // coprocessorinitializewithpython(sname, len);
    this->coProcInit = true;

    this->m_proc = vtkCPProcessor::New();
    this->m_proc->Initialize();
    this->m_link = vtkLiveInsituLink::New();
    this->m_link->SetInsituPort(22222);
    // this->m_link->SetInsituPort(11111);
    this->m_link->SetHostname("134.130.83.135");
    // this->m_link->SetHostname("localhost");
    // this->m_link->SetHostname("workstation00lnx");
    this->m_link->SetProcessType(vtkLiveInsituLink::INSITU);

    this->m_spxm =
        vtkSMProxyManager::GetProxyManager()->GetActiveSessionProxyManager();

    this->m_px = this->m_spxm->NewProxy("sources", "PVTrivialProducer");
    vtkSMProxy *px2 = this->m_spxm->NewProxy("sources", "PVTrivialProducer");

    this->m_spx = vtkSMSourceProxy::SafeDownCast(this->m_px);
    vtkSMSourceProxy *spx2 = vtkSMSourceProxy::SafeDownCast(px2);

    this->m_spxm->RegisterProxy("sources", "Mesh", this->m_spx);
    this->m_spxm->RegisterProxy("sources", "tplo", px2);

    vtkObjectBase *obase = this->m_spx->GetClientSideObject();
    vtkObjectBase *obase2 = spx2->GetClientSideObject();

    this->m_meshProd = vtkPVTrivialProducer::SafeDownCast(obase);
    this->m_tploProd = vtkPVTrivialProducer::SafeDownCast(obase2);

    this->m_meshProd->SetOutput(this->m_toplot, this->m_time);
    this->m_tploProd->SetOutput(this->m_tplo, this->m_time);

    this->m_meshProd->UpdateTimeStep(0.0);
    this->m_tploProd->UpdateTimeStep(0.0);

    this->m_link->Initialize(this->m_spxm);
  }

  this->Update();
#endif
}

/**
                finalizeCoProcessing closes the network connection to Paraview.
**/
inline void managementClass::finalizeCoProcessing() {
#ifdef USE_VTK
  if (this->coProcInit) {
    this->m_link->DropLiveInsituConnection();
  }
#endif
}

/**
                Update: Sends the current data to Paraview.
**/
inline void managementClass::Update() {
#ifdef USE_VTK
  if (this->coProcInit) {
    this->m_link->InsituUpdate(this->m_time, this->m_tstep);
    this->m_spx->UpdatePipeline(this->m_time);
    this->m_link->InsituPostProcess(this->m_time, this->m_tstep);
  }
#endif
}

inline void managementClass::checkTploArray(const char *name) {
#ifdef USE_VTK
  if (this->m_tplo->GetColumnByName(name) == 0) {
    vtkSmartPointer<vtkDoubleArray> toAdd =
        vtkSmartPointer<vtkDoubleArray>::New();
    toAdd->SetName(name);
    vtkIdType numrows = this->m_tplo->GetNumberOfRows();
    for (auto i = 0; i < numrows; ++i) {
      toAdd->InsertNextTuple1(0.0);
    }
    this->m_tplo->AddColumn(toAdd);
  }
#endif
}

inline void managementClass::addTploValue(const char *name, double &value) {
#ifdef USE_VTK
  this->checkTploArray(name);

  if (this->m_tplo->GetNumberOfRows() < this->m_row) {
    this->m_tplo->InsertNextBlankRow();
  }
  this->m_tplo->SetValueByName(this->m_row - 1, name, value);
  // this->checkTploSize();
#endif
}

inline void managementClass::checkTploSize() {
#ifdef USE_VTK
  vtkIdType ncols = this->m_tplo->GetNumberOfColumns();
  for (auto i = 0; i < ncols; ++i) {
    vtkDoubleArray *arr =
        vtkDoubleArray::SafeDownCast(this->m_tplo->GetColumn(i));
    while (arr->GetNumberOfTuples() <= this->m_row) {
      // std::cout << arr->GetNumberOfTuples();
      arr->InsertNextTuple1(0.0);
    }
  }
#endif
}
