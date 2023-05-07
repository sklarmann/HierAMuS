// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once


#include <plot/vtkplotClassBase.h>
#include <vector>
#include <types/MatrixTypes.h>

#include <memory>

class managementClass;

namespace HierAMuS {


class vtkPlotInterface : public vtkPlotInterfaceBase {
public:
  vtkPlotInterface();
  virtual ~vtkPlotInterface() = default;


  vtkUnstructuredGrid *getGrid() override;
  void interact() override;

  void timeUpdate(const prec &time) override;

  void initFileNames(std::string folder, std::string file) override;
  void toFile(std::string folder, std::string file) override;




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
  void initialize(std::string fileName) override;
  /**
   * @brief Starts CoProcessing with Paraview on localhost.
   *
   */
  void CoProcess() override;
  /**
   * @brief Finalizes CoProcessing.
   *
   */
  void finalizeCoProcessing() override;
  /**
   * @brief Switches to special routines in the Paraview management class to be an RVE.
   *
   */
  void isRVEfunc() override;
  /**
   * @brief Updates the time step in the Paraview management class
   *
   * @param[in] time Current time value.
   */
  void TimeUpdate(prec time) override;
  /**
   * @brief Writes the Paraview data to a file.
   *
   */
  void writeFile() override;
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
  void addPoint(indexType main, indexType part, indexType id,
                prec x, prec y, prec z) override;
  /**
   * @brief Adds a point to the Paraview Adapter.
   *
   * @param[in] main The number of the main block.
   * @param[in] part The number of the part block.
   * @param[in] id The id of the point to add.
   * @param[in] x,y,z coordinates of the point.
   */
  void addPoint(indexType main, indexType part, indexType id,
                Types::Vector3<prec> &coors) override;
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
  void addCell(indexType main, indexType part, indexType FeapCellNumber,
               indexType FeapSubCellNumber, std::vector <indexType> FeapPoints, indexType numpts,
               int vtkNumber) override;
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
  void addCompleteFeapPointArray(indexType main, indexType part,
                                 std::vector <prec> data, indexType num_comp,
                                 indexType numPoints, std::string name) override;
  void addCompleteFeapPointArray(indexType main, indexType part,
                                 std::vector <indexType> data, indexType num_comp,
                                 indexType numPoints, std::string name) override;
  /**
   * @brief Adds point data to all meshes of main.
   *
   * @param[in] main Main mesh number.
   * @param[in] data Data to add.
   * @param[in] num_comp Number of components per point.
   * @param[in] numPoints Number of points.
   * @param[in] name Name of the field.
   */
  template<typename T>
  void addCompleteFeapPointArrayToAll(indexType main, const std::vector <T> &data,
                                      indexType num_comp, indexType numPoints,
                                      const std::string &name);

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
  void setPointData(indexType main, indexType part, indexType feapId,
                    std::vector <indexType> data, indexType num_comp, std::string name) override;
  void setPointData(indexType main, indexType part, indexType feapId,
                    std::vector <prec> data, indexType num_comp, std::string name) override;

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
  void SumPointDataWeighted(indexType main, indexType part, std::vector <prec> &data,
                    indexType feapId, indexType num_comp, const std::string &name) override;
  void SumPointDataWeighted(indexType main, indexType part, std::vector <indexType> &data,
                    indexType feapId, indexType num_comp, const std::string &name) override;

  void SumPointData(indexType main, indexType part, std::vector <prec> &data,
                    indexType feapId, indexType num_comp, const std::string &name) override;
  void SumPointData(indexType main, indexType part, std::vector <indexType> &data,
                    indexType feapId, indexType num_comp, const std::string &name) override;

  std::string DisplacementName = "Displacements";
  std::string RotationName = "Rotations";
  std::string StressName = "Stresses";
  std::string StrainName = "Strains";
private:
  void pvdFileReader(std::string &pvdFile, std::vector <std::string> &FNames, std::vector <prec> &timesteps);
  void pvdFileWriter(std::string &pvdFile, std::vector <std::string> &FNames, std::vector <prec> &timesteps);

  // vtkSmartPointer<vtkPoints> points;
  // vtkSmartPointer<vtkCellArray> cells;
  // vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
  // vtkSmartPointer<vtkDataSetMapper> mapper;
  // vtkSmartPointer<vtkActor> actor;
  // vtkSmartPointer<vtkRenderer> renderer;
  // vtkSmartPointer<vtkRenderWindow> renderWindow;
  // vtkSmartPointer<vtkScalarBarActor> scalarBar;
  // vtkSmartPointer<vtkLookupTable> hueLut;
  indexType step;
  std::vector <prec> timesteps;
  std::shared_ptr<managementClass> paraViewManager;
//  managementClass paraViewManager;
  bool adapterInit;
};


}
