// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once



#include "Nodetypes.h"

#include <Eigen/Dense>
#include <map>
#include <vector>

#include <types/MatrixTypes.h>

#include <plot/vtkplotClass.h>
#include <plot/vtkplotClassBase.h>


#include "geometry/GeometryShape.h"

#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

#include "NormTypes.h"

#include "solver/HistoryDataNew/HistoryDataStructure.h"
#include "solver/HistoryDataNew/HistoryDataIterator.h"
#include "solver/HistoryDataNew/HistoryDataSetup.h"

#include "finiteElements/ElementTypes.h"

template <class bla> class vtkSmartPointer;

class vtkCell;

class managementClass;

namespace HierAMuS::Materials {
class Material;
}

namespace HierAMuS::Geometry {
class FacesData;
class VolumesData;
}

namespace HierAMuS::FiniteElement {


class GenericFiniteElement {
  using ptrCol = HierAMuS::PointerCollection;

public:
  GenericFiniteElement();
  virtual ~GenericFiniteElement();
  virtual auto getElementType() -> Elementtypes {
    return Elementtypes::Generic;
  };

  virtual void set_pointers(PointerCollection &pointers) = 0;


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

  void setId(indexType id) { this->m_id = id; };
  auto getId() -> indexType { return this->m_id; };
  void setMatrial(Materials::Material *in) { this->m_material = in; };
  virtual void setMaterialPerSubElement(std::vector<indexType> &materialNumbers){};
  auto getMaterial() -> HierAMuS::Materials::Material * {
    return this->m_material;
  };
  auto getMaterialFormulation(PointerCollection &pointers)
      -> std::shared_ptr<HierAMuS::Materials::GenericMaterialFormulation>;

  virtual auto getMaterialFormulation(PointerCollection &pointers,
                                      IntegrationPoint &ip)
      -> std::shared_ptr<HierAMuS::Materials::GenericMaterialFormulation>;
  auto getMaterialId() -> indexType;


  virtual void GenericSetDegreesOfFreedom(PointerCollection& pointers) = 0;
  virtual void GenericAdditionalOperations(PointerCollection& pointers) = 0;

  virtual auto getDofs(PointerCollection& pointers) -> std::vector<DegreeOfFreedom *>;

  virtual void GenericSetTangentResidual(
    PointerCollection& pointers,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) = 0;


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
      -> Geometry::VertexData &;
  virtual auto getEdge(ptrCol &pointers, indexType localNumber)
      -> Geometry::EdgesData &;
  virtual auto getFace(ptrCol &pointers, indexType localNumber)
      -> Geometry::FacesData * {
    throw std::runtime_error("Method getFace not implemented for Element");
    return nullptr;
  };
  virtual auto getVolume(ptrCol &pointers, indexType localNumber)
      -> Geometry::VolumesData * {
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



  
  virtual void getElementsLocalNodalReactions(
      ptrCol &pointers, std::map<indexType, std::vector<prec>> &vReacs);

  virtual void getSolution(ptrCol &pointers,
                           std::vector<DegreeOfFreedom *> &Dofs,
                           Types::VectorX<prec> &solution);
  auto getSolution(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs)
      -> Types::VectorX<prec>;
  auto getIncrementalSolution(ptrCol &pointers,
                              std::vector<DegreeOfFreedom *> &Dofs)
      -> Types::VectorX<prec>;;
  auto getNewtonSolution(ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs)
      -> Types::VectorX<prec>;;
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


  //L2Shapes
  virtual void getL2Dofs(ptrCol& pointers, std::vector<DegreeOfFreedom*>& Dofs,
      indexType meshID, indexType order) {
      throw std::runtime_error("getL2Dofs not implemented!");
  };
  virtual void setL2Shapes(ptrCol& pointers, indexType meshid,
      indexType order) {
      throw std::runtime_error("setL2Shapes not implemented!");
  };



  auto getNumberOfIntegrationPoints(PointerCollection& pointers) -> indexType;
  auto getElementHistoryDataStructure(PointerCollection& pointers) -> const HistoryDataStructure &;
  auto getMaterialHistoryDataStructure(PointerCollection& pointers) -> const HistoryDataStructure &;
  auto getHistoryDataSetUp(PointerCollection& pointers) -> HistoryDataSetup;
  auto getHistoryDataIterator(ptrCol &pointers) -> HistoryDataIterator;

  void updateRVEHistory(PointerCollection &pointers);

  // plot
  virtual void toParaviewAdapter(ptrCol &pointers, vtkPlotInterface &catalyst,
                         ParaviewSwitch ToDo) = 0;

  
  // Element data fields
  void request_element_data_field(PointerCollection &pointers,
                                  indexType fieldId,
                                  indexType rows, indexType cols);
  auto get_element_data_field(PointerCollection &pointers, 
                              indexType fieldId) -> Types::MatrixXX<prec> &;
  void set_element_data_field(PointerCollection &pointers, 
                              indexType fieldId, Types::MatrixXX<prec> &data);

protected:
  indexType m_id;
  HierAMuS::Materials::Material *m_material;
};
} // namespace HierAMuS::FiniteElement
