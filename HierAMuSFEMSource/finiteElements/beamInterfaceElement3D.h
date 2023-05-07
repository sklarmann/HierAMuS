// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "datatypes.h"
#include "equations/Nodetypes.h"
#include "geometry/BeamInterface3D.h"
#include "geometry/GeometryTypes.h"
#include "pointercollection/pointercollection.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <forwarddeclaration.h>

#include <finiteElements/GenericFiniteElement.h>
#include <types/MatrixTypes.h>

#include <Eigen/Dense>
#include <map>
#include <vector>
#include <set>

template <class bla> class vtkSmartPointer;

class vtkCell;

namespace HierAMuS::FiniteElement {

class beamInterfaceElement3D : public FiniteElement::GenericFiniteElement {
  using ptrCol = PointerCollection;

public:
  beamInterfaceElement3D();
  ~beamInterfaceElement3D() override;
  auto getElementType() -> Elementtypes override {
    return Elementtypes::beamInterfaceElement3D;
  };

  /**
   * @brief For each face of the cross-section, assigns the corresponding material number.  
   * 
   * @param materials 
   */
  void setMaterialPerSubElement(std::vector<indexType> &materials) override;

  
  /**
   * @brief Assigns the necessary degrees of freedom to the geometric elements.
   * 
   * @param pointers 
   * @param dispOrder 
   * @param meshIdDisp 
   * @param meshIdRot 
   */
  void setH1Shapes(PointerCollection &pointers, indexType dispOrder,
                   indexType meshIdDisp, indexType meshIdRot);

  /**
   * @brief Gets the H1 shapes associated with the corresponding degrees of freedom.
   * 
   * @param pointers 
   * @param dispOrder 
   * @param meshIdDisp 
   * @param jacobi 
   * @param integrationPoint 
   * @return Geometry::H1Shapes 
   */
  auto getH1Shapes(PointerCollection &pointers, indexType dispOrder,
                   indexType meshIdDisp, const Types::Matrix33<prec> &jacobi,
                   IntegrationPoint &integrationPoint) -> Geometry::H1Shapes;

  /**
   * @brief Get the Integration Points for the element.
   * 
   * @param pointers 
   * @return IntegrationPoints 
   */
  auto getIntegrationPoints(ptrCol &pointers) -> IntegrationPoints override;

  /**
   * @brief Computes the warping parameters for the element.
   *
   * The warping parameters are necessary for the approximation of the warping strains.
   * Shape functions are derived by assumed superposed displacements u_1, u_2, u_3.
   * There are several versions of the approximation:
   * - 1: Only special considerations are done for the normal strains.
   * - 2: 
   * 
   * @param pointers 
   */
  void computeWarpingShapes(PointerCollection &pointers);

  auto getMaterialFormulation(PointerCollection& pointers, IntegrationPoint &ip)
  -> std::shared_ptr<HierAMuS::Materials::GenericMaterialFormulation> override;

  struct warpingShapes3D {
    Types::VectorX<prec> omega1, omega2, omega3, localShapes;
    Types::Matrix3X<prec> omega1Deriv, omega2Deriv, omega3Deriv, localShapesDeriv;
  };
  auto getWarpingShapes(PointerCollection &pointers,
                        const Types::Matrix33<prec> &jacobi,
                        IntegrationPoint &integrationPoint)
      -> warpingShapes3D;

  auto getJacobian(PointerCollection &pointers,
                   IntegrationPoint &integrationPoint)
      -> Types::MatrixXX<prec> override;

  void setVerts(std::vector<indexType> &vertsIn) override;
  void setFaces(const std::vector<indexType> &facesIn) override;

  auto getH1Dofs(ptrCol &pointers, indexType meshIDdisp, indexType meshIDrot,
                 indexType order) -> std::vector<DegreeOfFreedom *>;

  auto getLocalCoordinate(PointerCollection &pointers,
                          IntegrationPoint &integrationPoint)
      -> Types::Vector3<prec>;
  auto getRotationR0() -> Types::Matrix33<prec>;


  struct localStressStrainInterface{
    indexType pos;
    indexType numShapes;
    Types::MatrixXX<prec> shapes;
  };
  auto getL2ShapesStressStrain(PointerCollection &pointers,
                               IntegrationPoint &ip, indexType order,
                               Types::MatrixXX<prec> &jacobi)
      -> localStressStrainInterface;

  auto getNumberOfInnerDofs(indexType order) -> indexType;

  
  void setWarpingType(indexType typ);
  auto getNumberOfInnerWarpingDofs() -> indexType;
  auto getWarpingMatrix(PointerCollection &pointers, IntegrationPoint &ip,
                        indexType order, Types::MatrixXX<prec> &jacobi)
      -> Types::Matrix6X<prec>;

  void print(PointerCollection &pointers);

  struct beamShapesSpecial{

  };
  auto getBeamShapesSpecial(PointerCollection &pointers, IntegrationPoint &ip,
                            indexType order, Types::MatrixXX<prec> &jacobi)
      -> beamShapesSpecial;

  auto getBMatrixWarpingLocal(PointerCollection &pointers, IntegrationPoint &ip,
                              indexType meshId, indexType order,
                              Types::MatrixXX<prec> &jacobi)
      -> Types::Matrix6X<prec>;
  auto getBMatrixWarpingLocal2(PointerCollection &pointers, IntegrationPoint &ip,
                              indexType meshId, indexType order,
                              Types::MatrixXX<prec> &jacobi)
      -> Types::Matrix6X<prec>;

  auto getWarpDofConstraint(PointerCollection& pointers) -> std::set<indexType>;
  void setWarpDofConstraintVertices(std::vector<indexType> &vertnumbers);


  // Paraview
  void geometryToParaview(PointerCollection &pointers,
                          vtkPlotInterface &paraviewAdapter, indexType mainMesh,
                          indexType subMesh) override;
  void computeWeightsParaview(PointerCollection &pointers,
                              vtkPlotInterface &paraviewAdapter,
                              indexType mainMesh, indexType subMesh) override;
  void H1SolutionToParaview(PointerCollection &pointers,
                            vtkPlotInterface &paraviewAdapter,
                            indexType mainMesh, indexType subMesh,
                            indexType meshId, indexType order,
                            std::string name) override;
  void projectDataToParaviewVertices(
      PointerCollection &pointers, vtkPlotInterface &paraviewAdapter,
      indexType mainMesh, indexType subMesh, indexType order,
      IntegrationPoint &IntegrationPt, Types::VectorX<prec> &data,
      indexType numberComponents, std::string name) override;



private:
  auto getElement(PointerCollection &pointers) -> Geometry::BeamInterface3D *;
  void createNodeShapeMapping(PointerCollection &pointers);
  void computeGeometry(PointerCollection &pointers);
  auto static getFaceIntegrationPoint(IntegrationPoint &integrationPoint)
      -> IntegrationPoint;
  void computeFunctionParameters(PointerCollection &pointers);

  auto getL2TransformationMatrix(PointerCollection &pointers,
                                 IntegrationPoint &ip) -> Types::Matrix66<prec>;

  auto getNumberOfInnerWarpingDofsV1() -> indexType;
  auto getWarpingMatrixV1(PointerCollection &pointers, IntegrationPoint &ip,
                          indexType order, Types::MatrixXX<prec> &jacobi)
      -> Types::Matrix6X<prec>;
  auto getNumberOfInnerWarpingDofsV2() -> indexType;
  auto getWarpingMatrixV2(PointerCollection &pointers, IntegrationPoint &ip,
                          indexType order, Types::MatrixXX<prec> &jacobi)
      -> Types::Matrix6X<prec>;
  auto getNumberOfInnerWarpingDofsV3() -> indexType;
  auto getWarpingMatrixV3(PointerCollection &pointers, IntegrationPoint &ip,
                          indexType order, Types::MatrixXX<prec> &jacobi)
      -> Types::Matrix6X<prec>;
  auto getNumberOfInnerWarpingDofsV4() -> indexType;
  auto getWarpingMatrixV4(PointerCollection &pointers, IntegrationPoint &ip,
                          indexType order, Types::MatrixXX<prec> &jacobi)
      -> Types::Matrix6X<prec>;

  void computeBCVerts(PointerCollection &pointers);


  // Warping shape data
  std::map<indexType, indexType> m_nodeShapeMapping;
  indexType m_H1MeshId; // Displacement nodes will be used to sort the warping shapes
  indexType m_MeshIdRot;
  indexType m_warpOrder;
  indexType m_numberOfWarpingShapes;
  Types::Matrix3X<prec> m_warpingCoefficients;
  Types::Matrix2X<prec> m_warpingCoefficientsUy, m_warpingCoefficientsUz;
  indexType m_warpingType;
  Types::Matrix3X<prec> m_surfToBeamShapesU;
  Types::Matrix3X<prec> m_surfToBeamShapesBeta; // Will contain the nodal values of the surface dofs to approximate a constant translation and rotation of the cross-section

  // Geometry data
  prec m_length;
  Types::Vector3<prec> m_projectedCoordinate;
  Types::Vector3<prec> m_A1;
  Types::Vector3<prec> m_A2;
  Types::Vector3<prec> m_A3;
  indexType m_shapeXiFace, m_shapeXiBeam; 

  // Description
  std::vector<indexType> m_faces;
  std::vector<indexType> m_materialPerFace;
  indexType m_beamVertex;

  indexType m_geoNumber;

  indexType m_vertMain, m_vertX2, m_vertX3;
  std::vector<indexType> m_warpEqIds;

  indexType minFaceVertNumber, maxFaceVertNumber;
};

} // namespace HierAMuS::FiniteElement
