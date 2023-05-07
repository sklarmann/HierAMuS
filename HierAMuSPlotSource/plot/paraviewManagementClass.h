// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#ifdef USE_VTK
#include <plot/catalystcxxinterface/management.h>
#endif


class paraviewManagementClass {
public:
  paraviewManagementClass();
  ~paraviewManagementClass();

  // Deprecated

  void addPoints(double *coor, int &numpoints, int &ndm);
#ifdef USE_VTK
  vtkUnstructuredGrid *getGrid(const int &main, const int &part);
#endif
  //

  void initialize(char *sname);
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
   
};