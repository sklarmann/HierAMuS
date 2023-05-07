// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>
#include <string>
#include <vector>

#include <vtkUnstructuredGrid.h>

#include <vtkAbstractArray.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>

static char pointIdName[12] = "Node Number";
static char cellIdName[15] = "Element Number";
static std::string weightName = "weights";

class gridHandler {
public:
  gridHandler();
  ~gridHandler();
  vtkUnstructuredGrid *getGrid();

  void addPoint(const int &FeapId, const double &x, const double &y,
                const double &z);

  template <typename T>
  void CompleteFeapPointArray(const T *data, const int &num_comp,
                              const int &numPoints, const char *name);
  template <typename T>
  void SetPointData(const T *data, const int &feapId, const int &num_comp,
                    const char *name);
  template <typename T>
  void SumPointDataWeighted(const T *data, const int &feapId, const int &num_comp,
                    const char *name);
  template <typename T>
  void SumPointData(const T *data, const int &feapId, const int &num_comp,
                    const char *name);
  template <typename T>
  void SetFieldData(const T *data, const int &num_comp, const char *name);

  void addCell(int FeapCellNumber, int FeapCellSubNumber, int *pointIds,
               int numpoints, int vtkNumber);

  template <typename T>
  void SetCellData(const T *data, const int &num_comp,
                   const int &feapCellNumber, const int &feapSubCellNumber,
                   const char *name);

  void setAllToZero();

  void normalizeFieldByWeight();

private:
  enum DataType { PointData, CellData };

  template <typename T>
  void createField(const char *name, const int &numComp, const DataType &type);

  bool hasField(const char *name, const DataType &type);

  vtkUnstructuredGrid *grid;

  std::map<int, vtkIdType>
      FeapToParvMapPoints; // Maps Feap Coor Id to Paraview Node Id
  std::map<int, std::vector<int>>
      FeapToParvMapCells; // Maps Feap Element Id to Paraview CellIds
  std::vector<std::string> arraysNormalizationNames;
};

inline gridHandler::gridHandler() {
  this->grid = vtkUnstructuredGrid::New();
  vtkPoints *points = vtkPoints::New();
  this->grid->SetPoints(points);
  vtkIntArray *feapIds = vtkIntArray::New();
  feapIds->SetName(pointIdName);
  feapIds->SetNumberOfComponents(1);
  vtkIntArray *feapCellIds = vtkIntArray::New();
  feapCellIds->SetName(cellIdName);
  feapCellIds->SetNumberOfComponents(2);
  this->grid->GetPointData()->AddArray(feapIds);
  this->grid->GetCellData()->AddArray(feapCellIds);
}

inline gridHandler::~gridHandler() {}

inline vtkUnstructuredGrid *gridHandler::getGrid() { return this->grid; }

template <typename T>
inline void
gridHandler::CompleteFeapPointArray(const T *data, const int &num_comp,
                                    const int &numPoints, const char *name) {

  if (!this->hasField(name, DataType::PointData)) {

    this->createField<T>(name, num_comp, DataType::PointData);
  }

  if (this->hasField(name, DataType::PointData)) {

    vtkDataArray *arr;
    int numnp = this->grid->GetNumberOfPoints();
    arr = this->grid->GetPointData()->GetArray(name);
    vtkDataArray *feapIds = vtkDataArray::SafeDownCast(
        this->grid->GetPointData()->GetArray(pointIdName));

    for (auto i = 0; i < numnp; ++i) {
      int feapId = feapIds->GetComponent(i, 0);
      for (auto j = 0; j < num_comp; ++j) {
        int pos = (feapId - 1) * num_comp;
        //this->grid->GetPointData()->GetArray(name)->SetComponent(i, j,
        //                                                         data[pos + j]);
        arr->SetComponent(i, j, data[pos]);
      }
    }
  }
}

template <typename T>
inline void gridHandler::SetPointData(const T *data, const int &feapId,
                                      const int &num_comp, const char *name) {
  if (!this->hasField(name, DataType::PointData)) {
    this->createField<T>(name, num_comp, DataType::PointData);
  }
  //
  auto search = this->FeapToParvMapPoints.find(feapId);
  if (search != this->FeapToParvMapPoints.end()) {
    vtkIdType vtkPointNumber = search->second;
    vtkDataArray *arr = this->grid->GetPointData()->GetArray(name);
    for (auto i = 0; i < num_comp; ++i) {
      arr->SetComponent(vtkPointNumber, i, data[i]);
    }
  }
}

template <typename T>
inline void gridHandler::SumPointDataWeighted(const T *data, const int &feapId,
                                      const int &num_comp, const char *name) {
  vtkDataArray *arr;
  if (!this->hasField(name, DataType::PointData)) {
    this->createField<T>(name, num_comp, DataType::PointData);
    arr = this->grid->GetPointData()->GetArray(name);
    arr->Fill(0.0);
    if (strcmp(name, weightName.c_str()) != 0) {
      this->arraysNormalizationNames.emplace_back(name);
    }
  } else {
    arr = this->grid->GetPointData()->GetArray(name);
  }

  //
  auto search = this->FeapToParvMapPoints.find(feapId);
  if (search != this->FeapToParvMapPoints.end()) {
    vtkIdType vtkPointNumber = search->second;
    T tempVal;
    for (auto i = 0; i < num_comp; ++i) {
      tempVal = arr->GetComponent(vtkPointNumber, i);
      tempVal += data[i];
      arr->SetComponent(vtkPointNumber, i, tempVal);
    }
  }
}

template <typename T>
inline void gridHandler::SumPointData(const T *data, const int &feapId,
                                      const int &num_comp, const char *name) {
  vtkDataArray *arr;
  if (!this->hasField(name, DataType::PointData)) {
    this->createField<T>(name, num_comp, DataType::PointData);
    arr = this->grid->GetPointData()->GetArray(name);
    arr->Fill(0.0);
  } else {
    arr = this->grid->GetPointData()->GetArray(name);
  }

  //
  auto search = this->FeapToParvMapPoints.find(feapId);
  if (search != this->FeapToParvMapPoints.end()) {
    vtkIdType vtkPointNumber = search->second;
    T tempVal;
    for (auto i = 0; i < num_comp; ++i) {
      tempVal = arr->GetComponent(vtkPointNumber, i);
      tempVal += data[i];
      arr->SetComponent(vtkPointNumber, i, tempVal);
    }
  }
}

template <typename T>
inline void gridHandler::SetFieldData(const T *data, const int &num_comp,
                                      const char *name) {
  vtkDataArray *arr;
  if (!this->grid->GetFieldData()->GetArray(name)) {
    if (std::is_same<T, int>()) {
      vtkIntArray *toadd = vtkIntArray::New();
      toadd->SetName(name);
      toadd->SetNumberOfComponents(num_comp);
      for (auto i = 0; i < num_comp; ++i)
        toadd->InsertNextValue(0);
      this->grid->GetFieldData()->AddArray(toadd);
    } else if (std::is_same<T, double>()) {
      vtkDoubleArray *toadd = vtkDoubleArray::New();

      toadd->SetName(name);
      toadd->SetNumberOfComponents(num_comp);
      for (auto i = 0; i < num_comp; ++i)
        toadd->InsertNextValue(0);
      this->grid->GetFieldData()->AddArray(toadd);
    }
  }
  arr = this->grid->GetFieldData()->GetArray(name);
  for (auto i = 0; i < num_comp; ++i) {
    arr->SetComponent(0, i, data[i]);
  }
}

template <typename T>
inline void gridHandler::createField(const char *name, const int &numComp,
                                     const DataType &type) {

  bool present = true;
  if (type == DataType::PointData) {
    present = this->hasField(name, DataType::PointData);
  } else {
    present = this->hasField(name, DataType::CellData);
  }

  if (!present) {
    vtkAbstractArray *toAdd;
    if (std::is_same<T, int>()) {
      toAdd = vtkIntArray::New();
    } else if (std::is_same<T, double>()) {
      toAdd = vtkDoubleArray::New();
    }
    toAdd->SetName(name);
    toAdd->SetNumberOfComponents(numComp);
    if (type == DataType::PointData) {
      toAdd->SetNumberOfTuples(this->grid->GetNumberOfPoints());
      this->grid->GetPointData()->AddArray(toAdd);
    } else {
      toAdd->SetNumberOfTuples(this->grid->GetNumberOfCells());
      this->grid->GetCellData()->AddArray(toAdd);
    }
  }
}

inline void gridHandler::addPoint(const int &FeapId, const double &x,
                                  const double &y, const double &z) {
  auto search = this->FeapToParvMapPoints.find(FeapId);
  if (search != this->FeapToParvMapPoints.end()) {

  } else {
    this->grid->GetPoints()->InsertNextPoint(x, y, z);
    vtkIdType parvnum = this->grid->GetPoints()->GetNumberOfPoints();
    this->FeapToParvMapPoints[FeapId] = parvnum - 1;
    this->grid->GetPointData()
        ->GetArray(pointIdName)
        ->InsertComponent(parvnum - 1, 0, FeapId);
  }
}

inline void gridHandler::addCell(int FeapCellNumber, int FeapCellSubNumber,
                                 int *pointIds, int numpoints, int vtkNumber) {

  auto vecSize = this->FeapToParvMapCells[FeapCellNumber].size();

  while (vecSize < FeapCellSubNumber) {
    this->FeapToParvMapCells[FeapCellNumber].emplace_back(-1);
    ++vecSize;
  }

  if (this->FeapToParvMapCells[FeapCellNumber][FeapCellSubNumber - 1] == -1) {

    vtkIdList *list = vtkIdList::New();
    for (auto i = 0; i < numpoints; ++i) {
      list->InsertNextId(this->FeapToParvMapPoints[pointIds[i]]);
    }
    this->grid->InsertNextCell(vtkNumber, list);
    int numcells = this->grid->GetNumberOfCells();
    this->grid->GetCellData()
        ->GetArray(cellIdName)
        ->InsertComponent(numcells - 1, 0, FeapCellNumber);
    this->grid->GetCellData()
        ->GetArray(cellIdName)
        ->InsertComponent(numcells - 1, 1, FeapCellSubNumber);
    this->FeapToParvMapCells[FeapCellNumber][FeapCellSubNumber - 1] =
        numcells - 1;
  }
}

inline bool gridHandler::hasField(const char *name, const DataType &type) {

  switch (type) {
  case DataType::PointData:
    if (this->grid->GetPointData()->GetArray(name)) {
      return true;
    }
    break;
  case DataType::CellData:
    if (this->grid->GetCellData()->GetArray(name)) {
      return true;
    }
    break;
  default:
    break;
  }

  return false;
}

inline void gridHandler::setAllToZero() {
  int numArrs = this->grid->GetPointData()->GetNumberOfArrays();
  vtkDataArray *arr;
  for (auto i = 0; i < numArrs; ++i) {
    arr = this->grid->GetPointData()->GetArray(i);
    std::string arrName(arr->GetName());
    std::string arrNameComp(pointIdName);
    if (arrName != arrNameComp) {
      arr->Fill(0.0);
    }
  }
}

template <typename T>
inline void gridHandler::SetCellData(const T *data, const int &num_comp,
                                     const int &feapCellNumber,
                                     const int &feapSubCellNumber,
                                     const char *name) {
  this->createField<T>(name, num_comp, DataType::CellData);

  auto search = this->FeapToParvMapCells.find(feapCellNumber);
  if (search != this->FeapToParvMapCells.end()) {
    if (feapSubCellNumber <= search->second.size()) {
      vtkIdType cellNum = search->second[feapSubCellNumber - 1];
      vtkDataArray *arr = this->grid->GetCellData()->GetArray(name);
      for (auto i = 0; i < num_comp; ++i) {
        arr->SetComponent(cellNum, i, data[i]);
      }
    }
  }
}

inline void gridHandler::normalizeFieldByWeight() {

  vtkDataArray *wArr = this->grid->GetPointData()->GetArray(weightName.c_str());
  if (wArr != NULL) {
    vtkIdType numComp;
    vtkIdType numTuples = wArr->GetNumberOfTuples();
    for (auto i = 0; i < this->arraysNormalizationNames.size(); ++i) {
      vtkDataArray *toWeight = this->grid->GetPointData()->GetArray(
          this->arraysNormalizationNames[i].c_str());
      numComp = toWeight->GetNumberOfComponents();
      for (vtkIdType j = 0; j < numTuples; ++j) {
        for (vtkIdType k = 0; k < numComp; ++k) {
          toWeight->SetComponent(
              j, k, toWeight->GetComponent(j, k) / wArr->GetComponent(j, 0));
          // std::cout << "\n" << j << " " << k << std::endl;
          // std::cout << toWeight->GetComponent(j,k) << " " <<
          // wArr->GetComponent(j,k) << " " <<
          // toWeight->GetComponent(j,k)/wArr->GetComponent(j,k) << " " <<
          // toWeight << " " << wArr << std::endl;
        }
      }
    }
  }
}
