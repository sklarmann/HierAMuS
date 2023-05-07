// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <forwarddeclaration.h>
#include <datatypes.h>
#include <memory>
#include "types/MatrixTypes.h"

class vtkUnstructuredGrid;
class vtkRenderer;

namespace HierAMuS {
enum class ParaviewSwitch {
  Mesh, Solution, ProjectedValues, Eigenvector, Weights
};

class paraviewNames {
public:
  static std::string weightName() { return "weights"; };
  static std::string stressName() { return "$\\sigma$"; };
  static std::string stressName9() { return "$\\sigma^9$"; };
  static std::string strainName() { return "$\\varepsilon$"; };
  static std::string strainName9() { return "$\\varepsilon^9$"; };
  static std::string DisplacementName() { return "u"; };
  static std::string TemperatureName() { return "Temperature"; };
  static std::string TemperatureGradientName() { return "TemperatureGradien"; };
  static std::string HeatFluxName() { return "HeatFlux"; };
  static std::string RotationName() { return "$\\omega$"; };
  static std::string ReactionsTransName() { return "ReactionsTrans"; };
  static std::string ReactionsRotName() { return "ReactionsRot"; };
};


class vtkPlotInterfaceBase {
public:
  typedef std::shared_ptr<vtkPlotInterfaceBase> sharedPtr;


  virtual vtkUnstructuredGrid *getGrid() { return nullptr; };
  virtual void interact(){};

  virtual void timeUpdate(const prec &time){};

  virtual void initFileNames(std::string folder, std::string file){};
  virtual void toFile(std::string folder, std::string file){};


  // Interface to Management Class
  /**
   * @brief Initialiazes the Paraview Management Class.
   * Takes the filename of the file including the directory path.
   * Extracts the filename and path. First letter of filename is replaced by R.
   * Append the folder parvout to the path.
   * Ouputfiles for Paraview are written into this
   *
   * @param[in] fileName File name of the FE program input file including paths.
   */
  virtual void initialize(std::string fileName){};
  /**
   * @brief Starts CoProcessing with Paraview on localhost.
   *
   */
  virtual void CoProcess(){};
  /**
   * @brief Finalizes CoProcessing.
   *
   */
  virtual void finalizeCoProcessing(){};
  /**
   * @brief Switches to special routines in the Paraview management class to be
   * an RVE.
   *
   */
  virtual void isRVEfunc(){};
  /**
   * @brief Updates the time step in the Paraview management class
   *
   * @param[in] time Current time value.
   */
  virtual void TimeUpdate(prec time){};
  /**
   * @brief Writes the Paraview data to a file.
   *
   */
  virtual void writeFile(){};
  /**
   * @brief Adds a point to the Paraview Adapter.
   *
   * @param[in] main The number of the main block.
   * @param[in] part The number of the part block.
   * @param[in] id The id of the point to add.
   * @param[in] x x-coordinate of the point.
   * @param[in] y y-coordinate of the point.
   * @param[in] z z-coordinate of the point.
   */
  virtual void addPoint(indexType main, indexType part, indexType id, prec x,
                        prec y, prec z){};
  /**
   * @brief Adds a point to the Paraview Adapter.
   *
   * @param[in] main The number of the main block.
   * @param[in] part The number of the part block.
   * @param[in] id The id of the point to add.
   * @param[in] x,y,z coordinates of the point.
   */
  virtual void addPoint(indexType main, indexType part, indexType id,
                        Types::Vector3<prec> &coors){};
  /**
   * @brief Adds a cell to the Paraview adapter
   *
   * @param[in] main Main mesh number.
   * @param[in] part Sub mesh number.
   * @param[in] FeapCellNumber Main element number.
   * @param[in] FeapSubCellNumber Sub element number.
   * @param[in] FeapPoints Point ids the cell is associated with.
   * @param[in] numpts Number of points the cell consists of.
   * @param[in] vtkNumber Vtk number of the cell type.
   */
  virtual void addCell(indexType main, indexType part, indexType FeapCellNumber,
                       indexType FeapSubCellNumber,
                       std::vector<indexType> FeapPoints, indexType numpts,
                       int vtkNumber){};
  /**
   * @brief Adds a point data field to the specific mesh.
   * Data are arranged in a vector of the length num_comp*numPoints.
   *
   * @param[in] main Main mesh number.
   * @param[in] part Part mesh number.
   * @param[in] data Vector of data to add.
   * @param[in] num_comp Number of components per point.
   * @param[in] numPoints Number of points.
   * @param[in] name Field name of the data.
   */
  virtual void addCompleteFeapPointArray(indexType main, indexType part,
                                         std::vector<prec> data,
                                         indexType num_comp,
                                         indexType numPoints,
                                         std::string name){};
  virtual void addCompleteFeapPointArray(indexType main, indexType part,
                                         std::vector<indexType> data,
                                         indexType num_comp,
                                         indexType numPoints,
                                         std::string name){};
  /**
   * @brief Adds point data to all meshes of main.
   *
   * @param[in] main Main mesh number.
   * @param[in] data Data to add.
   * @param[in] num_comp Number of components per point.
   * @param[in] numPoints Number of points.
   * @param[in] name Name of the field.
   */
  template <typename T>
  void addCompleteFeapPointArrayToAll(indexType main,
                                      const std::vector<T> &data,
                                      indexType num_comp, indexType numPoints,
                                      const std::string &name){};

  /**
   * @brief Sets data of a specific point.
   *
   * @param[in] main Main mesh number.
   * @param[in] part Part mesh number.
   * @param[in] feapId Point id.
   * @param[in] data Data to set.
   * @param[in] num_comp Number of components in data.
   * @param[in] name Name of the field.
   */
  virtual void setPointData(indexType main, indexType part, indexType feapId,
                            std::vector<indexType> data, indexType num_comp,
                            std::string name){};
  virtual void setPointData(indexType main, indexType part, indexType feapId,
                            std::vector<prec> data, indexType num_comp,
                            std::string name){};

  /**
   * @brief Sums up given data at the specific point.
   *
   * @param[in] main Main mesh number.
   * @param[in] part Part mesh number.
   * @param[in] data Data to sum up
   * @param[in] feapId Point id.
   * @param[in] num_comp Number of components in data.
   * @param[in] name Name of the field.
   */
  virtual void SumPointDataWeighted(indexType main, indexType part,
                                    std::vector<prec> &data, indexType feapId,
                                    indexType num_comp,
                                    const std::string &name){};
  virtual void SumPointDataWeighted(indexType main, indexType part,
                                    std::vector<indexType> &data,
                                    indexType feapId, indexType num_comp,
                                    const std::string &name){};

  virtual void SumPointData(indexType main, indexType part,
                            std::vector<prec> &data, indexType feapId,
                            indexType num_comp, const std::string &name){};
  virtual void SumPointData(indexType main, indexType part,
                            std::vector<indexType> &data, indexType feapId,
                            indexType num_comp, const std::string &name){};

};

}

