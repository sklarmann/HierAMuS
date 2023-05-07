// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <equations/DofStatus.h>
#include <finiteElements/ElementTypes.h>
#include <forwarddeclaration.h>
#include <geometry/GeometryData.h>
#include <materials/Material.h>

#include <Eigen/Dense>
#include <map>
#include <vector>

#include <types/MatrixTypes.h>

#include <plot/vtkplotClass.h>
#include <plot/vtkplotClassBase.h>

#include "equations/GenericNodes.h"
#include "geometry/Base.h"
#include "pointercollection/pointercollection.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

#include "NormTypes.h"

#include "solver/HistoryDataNew/HistoryDataStructure.h"
#include "solver/HistoryDataNew/HistoryDataIterator.h"
#include "solver/HistoryDataNew/HistoryDataSetup.h"

template <class bla> class vtkSmartPointer;

class vtkCell;

class managementClass;

namespace HierAMuS::FiniteElement {

class GenericFiniteElement {
  using ptrCol = HierAMuS::PointerCollection;

public:
  GenericFiniteElement();
  virtual ~GenericFiniteElement();
  virtual auto getElementType() -> Elementtypes {
    return Elementtypes::Generic;
  };


  // Local Basis
  virtual auto getA1Vector(ptrCol &pointers,
                           IntegrationPoint &integration_point)
      -> Types::Vector3<prec> {
    throw std::runtime_error(
        "Get A1 vector for the specific element not implemented!");
    return {};
  }

  virtual void setAllNodeBoundaryConditionMeshId(ptrCol &pointers,
                                                 indexType meshId,
                                                 indexType dof);

  virtual auto getIntegrationPoints(ptrCol &pointers) -> IntegrationPoints;

  void setId(indexType id) { this->id = id; };
  auto getId() -> indexType { return this->id; };
  void setMatrial(Materials::Material *in) { this->material = in; };
  virtual void setMaterialPerSubElement(std::vector<indexType> &materialNumbers){};
  auto getMaterial() -> HierAMuS::Materials::Material * {
    return this->material;
  };
  auto getMaterialFormulation(PointerCollection& pointers)
  -> std::shared_ptr<HierAMuS::Materials::GenericMaterialFormulation> {
    return this->material->getMaterialFormulation(pointers);
  };

  virtual auto getMaterialFormulation(PointerCollection& pointers, IntegrationPoint &ip)
  -> std::shared_ptr<HierAMuS::Materials::GenericMaterialFormulation> {
    return this->material->getMaterialFormulation(pointers);
  };
  auto getMaterialId() -> indexType { return this->material->getNumber(); }

  void insertStiffnessResidual(
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, 1, Eigen::Dynamic> &residual,
      std::vector<indexType> &eqIds, std::vector<dofStatus> &eqStatus);

  void GenericSetDegreesOfFreedom(PointerCollection& pointers);
  void GenericAdditionalOperations(PointerCollection& pointers);

  auto getDofs(PointerCollection& pointers) -> std::vector<DegreeOfFreedom *>;

  void GenericSetTangentResidual(
    PointerCollection& pointers,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);


  auto computeNorm(ptrCol &pointers, NormTypes type) -> prec;
  void
  GenericSetMass(PointerCollection& pointers,
                 Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
                 Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  // TODO throw exception
  virtual void setVerts(std::vector<indexType> &vertsIn){};
  virtual void setEdges(std::vector<indexType> &edgesIn){};
  virtual void setFace(indexType face){};
  virtual void setFaces(const std::vector<indexType> &facesIn){};
  virtual void setVolume(indexType face){};
  virtual void setSpecial(std::vector<indexType> &specialIn){};
  virtual void setSpecial(indexType specialIn){};
  virtual auto getType() -> Elementtypes { return Elementtypes::Generic; };
  virtual auto getVertex(ptrCol &pointers, indexType localNumber)
      -> Geometry::Vertex &;
  virtual auto getEdge(ptrCol &pointers, indexType localNumber)
      -> Geometry::Edges &;
  virtual auto getFace(ptrCol &pointers, indexType localNumber)
      -> Geometry::Faces * {
    throw std::runtime_error("Method getFace not implemented for Element");
    return nullptr;
  };
  virtual auto getVolume(ptrCol &pointers, indexType localNumber)
      -> Geometry::Volumes * {
    throw std::runtime_error("Method getVolume not implemented for Element");
    return nullptr;
  };
  virtual auto getVertexId(ptrCol &pointers, indexType num) -> indexType {
    return 0;
  };
  virtual void setStiffnessMatrix(){};
  virtual auto getVertexIds(PointerCollection& pointers) -> std::vector<indexType> {
    std::vector<indexType> u;
    u.resize(0);
    return u;
  };
  virtual void getNumberOfNodes(ptrCol &pointers, indexType &numEdge1,
                                indexType &numEdge2, indexType &numEdge3,
                                indexType &numFace, indexType meshId){};
  virtual void
  getVertexFunctions(Eigen::Matrix<prec, 1, Eigen::Dynamic> &shape,
                     Eigen::Matrix<prec, 2, Eigen::Dynamic> &shapeDerivative,
                     prec xsi, prec eta){};

  virtual void getNodesMeshId(ptrCol &pointers,
                              std::vector<GenericNodes *> &Nodes,
                              indexType meshId){};

  virtual void
  getAffineCoordinates(Eigen::Matrix<prec, 1, Eigen::Dynamic> &coor,
                       Eigen::Matrix<prec, 2, Eigen::Dynamic> &coorDerivative,
                       prec xsi, prec eta){};

  virtual void getJacobian(ptrCol &pointers, Eigen::Matrix<prec, 2, 2> &jacobi,
                           Eigen::Matrix<prec, 2, Eigen::Dynamic> &coorDeriv){};

  
  virtual void getElementsLocalNodalReactions(
      ptrCol &pointers, std::map<indexType, std::vector<prec>> &vReacs);

  virtual void getSolution(ptrCol &pointers,
                           std::vector<DegreeOfFreedom *> &Dofs,
                           Types::VectorX<prec> &solution);
  virtual auto getSolution(ptrCol &pointers,
                           std::vector<DegreeOfFreedom *> &Dofs)
      -> Types::VectorX<prec>;
  virtual void getVelocity(ptrCol &pointers,
                           std::vector<DegreeOfFreedom *> &Dofs,
                           Types::VectorX<prec> &solution);
  virtual void getAcceleration(ptrCol &pointers,
                               std::vector<DegreeOfFreedom *> &Dofs,
                               Types::VectorX<prec> &solution);

  // Intergration Points
  virtual void getGaussPoints(indexType number, std::vector<prec> &weight,
                              std::vector<prec> &xsi) {
    throw std::runtime_error("getGaussPoints 1D not implemented!");
  };
  virtual void getGaussPoints(indexType number, std::vector<prec> &weight,
                              std::vector<prec> &xsi, std::vector<prec> &eta) {
    throw std::runtime_error("getGaussPoints 2D not implemented!");
  };
  virtual void getGaussPoints(indexType number, std::vector<prec> &weight,
                              std::vector<prec> &xsi, std::vector<prec> &eta,
                              std::vector<prec> &zeta) {
    throw std::runtime_error("getGaussPoints 3D not implemented!");
  };

  // Structure
  /**
   * @brief getNumberOfVertices returns the number of vertices.
   * @return indexType, the number of vertices.
   */
  virtual auto getNumberOfVertices(PointerCollection& pointers) -> indexType {
    throw std::runtime_error("getNumberOfVertices not implemented!");
  };
  /**
   * @brief getNumberOfEdges returns the number of edges.
   * @return indexType, the number of edges.
   */
  virtual auto getNumberOfEdges(PointerCollection& pointers) -> indexType {
    throw std::runtime_error("getNumberOfEdges not implemented!");
  };

  /**
   * @brief getNumberOfFace returns the number of faces.
   * @return indexType, the number of faces.
   */
  virtual auto getNumberOfFaces(PointerCollection& pointers) -> indexType {
    throw std::runtime_error("getNumberOfFaces not implemented!");
  };

  // Geometric mapping
  virtual auto getJacobian(ptrCol &pointers, IntegrationPoint &IntegrationPt)
      -> Types::MatrixXX<prec> {
    std::runtime_error("getJacobian not implemented for this element!");
    return {};
  }
  virtual void getJacobian(ptrCol &pointers, prec &jacobi, prec xsi) {
    std::runtime_error("getJacobian 1D not implemented for this element!");
  };
  virtual void getJacobian(ptrCol &pointers, Types::Matrix22<prec> &jacobi,
                           prec xsi, prec eta) {
    std::runtime_error("getJacobian 2D not implemented for this element!");
  };
  virtual void getJacobian(ptrCol &pointers, Types::Matrix33<prec> &jacobi,
                           prec xsi, prec eta, prec zeta) {
    std::runtime_error("getJacobian 3D not implemented for this element!");
  };

  // H1 Shape Functions
  virtual void setH1Shapes(ptrCol &pointers, indexType meshid,
                           indexType order) {
    throw std::runtime_error("setH1Shapes not implemented!");
  };
  virtual void getH1Dofs(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs,
                         indexType meshID, indexType order) {
    throw std::runtime_error("getH1Dofs not implemented!");
  };
  virtual auto getH1Nodes(ptrCol &pointers, indexType meshID, indexType order)
      -> std::vector<GenericNodes *> {
    throw std::runtime_error("getH1Dofs not implemented!");
    return {};
  };
  virtual auto getH1Shapes(ptrCol &pointers, indexType order,
                           Types::MatrixXX<prec> &jacobi,
                           IntegrationPoint &IntegrationPt)
      -> Geometry::H1Shapes {
    throw std::runtime_error("getH1Shapes not implemented!");
    return {};
  };
  virtual void getH1Shapes(ptrCol &pointers, indexType order,
                           Types::MatrixXX<prec> &jacobi,
                           Types::VectorX<prec> &shape,
                           Types::MatrixXX<prec> &shapeDerivative,
                           IntegrationPoint &IntegrationPt) {
    throw std::runtime_error("getH1Shapes not implemented!");
  };
  virtual void getH1Shapes(ptrCol &pointers, indexType order, prec jacobi,
                           Types::VectorX<prec> &shape,
                           Types::VectorX<prec> &shapeDerivative, prec xsi) {
    throw std::runtime_error("getH1Shapes 1D not implemented!");
  };
  virtual void getH1Shapes(ptrCol &pointers, indexType order,
                           Types::Matrix22<prec> &jacobi,
                           Types::VectorX<prec> &shape,
                           Types::Matrix2X<prec> &shapeDerivative, prec xsi,
                           prec eta) {
    throw std::runtime_error("getH1Shapes 2D not implemented!");
  };
  virtual void getH1Shapes(ptrCol &pointers, indexType order,
                           Types::Matrix33<prec> &jacobi,
                           Types::VectorX<prec> &shape,
                           Types::Matrix3X<prec> &shapeDerivative, prec xsi,
                           prec eta, prec zeta) {
    throw std::runtime_error("getH1Shapes 3D not implemented!");
  };

  // HCurl Shapes

  // Special Plate shapes
  virtual void setSpecialPlateShapes(PointerCollection &pointers,
                                     indexType meshid, indexType order){};
  virtual auto getSpecialPlateDofs(PointerCollection &pointers,
                                   indexType meshID, indexType order)
      -> std::vector<DegreeOfFreedom *> {
    return {};
  };
  virtual auto getSpecialPlateShapes(PointerCollection &pointers,
                                     indexType order,
                                     Types::MatrixXX<prec> &jacobi,
                                     IntegrationPoint &intPoint)
      -> Geometry::SpecialPlateShapes {
    return {};
  };

  // HDiv Shapes
  virtual void setHDivShapes(ptrCol &pointers, indexType meshid,
                             indexType order, NodeTypes type) {
    throw std::runtime_error("setHDivShapes not implemented!");
  };
  virtual void getHDivDofs(ptrCol &pointers,
                           std::vector<DegreeOfFreedom *> &Dofs,
                           indexType meshID, indexType order) {
    throw std::runtime_error("getHDivDofs not implemented!");
  };
  virtual void getHDivShapes(ptrCol &pointers, indexType order,
                             Types::Matrix22<prec> &jacobi,
                             Types::Matrix2X<prec> &shape,
                             Types::VectorX<prec> &dshape, prec xi,
                             prec eta) {
    throw std::runtime_error("getHDivShapes not implemented!");
  };

  virtual auto getHDivShapes(PointerCollection &pointers, indexType order,
                             Types::MatrixXX<prec> &jacobi,
                             IntegrationPoint &IntegrationPt)
      -> Geometry::HDivShapes {
    throw std::runtime_error("getHDivShapes not implemented!");
    return {};
  };

  //L2Shapes
  virtual void getL2Dofs(ptrCol& pointers, std::vector<DegreeOfFreedom*>& Dofs,
      indexType meshID, indexType order) {
      throw std::runtime_error("getL2Dofs not implemented!");
  };
  virtual void setL2Shapes(ptrCol& pointers, indexType meshid,
      indexType order) {
      throw std::runtime_error("setL2Shapes not implemented!");
  };
  virtual auto getL2Shapes(ptrCol& pointers, indexType order,
      Types::MatrixXX<prec>& jacobi,
      IntegrationPoint& IntegrationPt)
      -> Geometry::L2Shapes {
      throw std::runtime_error("getL2Shapes not implemented!");
      return {};
  }


  auto getNumberOfIntegrationPoints(PointerCollection& pointers) -> indexType;
  auto getElementHistoryDataStructure(PointerCollection& pointers) -> const HistoryDataStructure &;
  auto getMaterialHistoryDataStructure(PointerCollection& pointers) -> const HistoryDataStructure &;
  virtual auto getHistoryDataSetUp(PointerCollection& pointers) -> HistoryDataSetup;
  virtual auto getHistoryDataIterator(ptrCol &pointers) -> HistoryDataIterator;

  void updateRVEHistory(PointerCollection &pointers);

  // plot
  void toParaviewAdapter(ptrCol &pointers, vtkPlotInterface &catalyst,
                         ParaviewSwitch ToDo);

  virtual void geometryToParaview(PointerCollection &pointers,
                                  vtkPlotInterface &paraviewAdapter,
                                  indexType mainMesh, indexType subMesh){};
  virtual void computeWeightsParaview(PointerCollection &pointers,
                                      vtkPlotInterface &paraviewAdapter,
                                      indexType mainMesh, indexType subMesh){};
  virtual void H1SolutionToParaview(PointerCollection &pointers,
                                    vtkPlotInterface &paraviewAdapter,
                                    indexType mainMesh, indexType subMesh,
                                    indexType meshId, indexType order,
                                    std::string name){};
  virtual void H1DataToParaview(PointerCollection &pointers,
                                vtkPlotInterface &paraviewAdapter,
                                indexType mainMesh, indexType subMesh,
                                Types::VectorX<prec> &Data,
                                indexType numberComponents, indexType order,
                                std::string &name){};
  virtual void projectDataToParaviewVertices(
      PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
      indexType mainMesh, indexType subMesh, indexType order,
      IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
      indexType numberOfComponents, std::string name){};

  // Deprecated
  virtual void getVtkCell(ptrCol &pointers, vtkSmartPointer<vtkCell> &cell){};

  

  virtual void projectOnVertsParaview(ptrCol &pointers,
                                      vtkPlotInterface &catalyst,
                                      Types::VectorX<prec> &values, prec &xsi,
                                      prec &eta, prec &weight,
                                      std::string name){};

  virtual void projectOnVertsParaview(ptrCol &pointers,
                                      vtkPlotInterface &catalyst,
                                      Types::VectorX<prec> &values, prec &xsi,
                                      prec &eta, prec &zeta, prec &weight,
                                      std::string name){};

  virtual void setParaviewCellData(ptrCol &pointers,
                                   vtkPlotInterface &catalyst);
  

protected:
  indexType id;
  indexType numberOfIntegrationPoints;

  HierAMuS::Materials::Material *material;
};
} // namespace HierAMuS::FiniteElement
