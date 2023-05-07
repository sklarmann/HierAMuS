// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <forwarddeclaration.h>

#include <Eigen/Dense>
#include <elementFormulations/GenericElementFormulation.h>
#include <types/MatrixTypes.h>

namespace HierAMuS::Elementformulations {


class EL104_TimoshenkoPrism
    : public GenericElementFormulation {
public:
  explicit EL104_TimoshenkoPrism(PointerCollection *ptrCol);
  ~EL104_TimoshenkoPrism() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) override;
  auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> std::vector<DegreeOfFreedom *> override;
  void setTangentResidual(
    PointerCollection& pointers,
    FiniteElement::GenericFiniteElement *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;

  void toParaviewAdaper(PointerCollection &pointers,
                                FiniteElement::GenericFiniteElement *elem,
                                vtkPlotInterface &paraviewAdapter,
                                ParaviewSwitch control) override;

  auto getHistoryDataStructure() -> const HistoryDataStructure & override;
  auto getNumberOfIntergrationPoints(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> indexType override;

private:
  void toParaviewPrism(PointerCollection &pointers,
                       FiniteElement::LinearPrism *elem,
                       vtkPlotInterface &paraviewAdapter,
                       ParaviewSwitch control);

  Types::Matrix66<prec> getMaterialTangent();
  Types::Matrix66<prec> getnonMaterialTangent();
  Types::Matrix6X<prec> getLocalStrainStressInterpolation();
  Types::Matrix3X<prec> Hmatrix(const Types::VectorX <prec>& solution, indexType numberBeamNodes);
  Types::Matrix3X<prec> Rmatrix(const Types::VectorX <prec>& solution, indexType numberBeamNodes);
  Types::Vector3<prec> xo_prime_matrix(const Types::VectorX <prec>& shapeDerivative,
                                       const Types::VectorX <prec>& solution,
                                       const Types::Matrix3X<prec>& coors,
                                       indexType numberBeamNodes);

  Types::Vector3<prec> Xo_prime_matrix(const Types::VectorX <prec>& shapeDerivative,
                                       const Types::Matrix3X<prec>& coors,
                                       indexType numberBeamNodes);

  Types::Matrix6X <prec> getBMatrix(const Types::VectorX <prec>& shapeDerivative,
                                    const Types::VectorX<prec>& shape,
                                    const Types::Matrix33<prec>& R0,
                                    indexType numDofs, indexType numberBeamNodes);

  Types::Matrix6X <prec> getNonBMatrix(const Types::VectorX <prec>& shapeDerivative,
                                       const Types::VectorX<prec>& shape,
                                       const Types::Matrix33<prec>& R01,
                                       const Types::Matrix33<prec>& R02,
                                       const Types::Matrix3X<prec>& coors,
                                       const Types::VectorX <prec>& solution,
                                       indexType numDofs, indexType numberBeamNodes);

  Types::Vector6 <prec> getstrain_Vector(const Types::VectorX <prec>& shapeDerivative,
                                         const Types::VectorX<prec>& shape,
                                         const Types::Matrix33<prec>& R01,
                                         const Types::Matrix33<prec>& R02,
                                         const Types::Matrix3X<prec>& coors,
                                         const Types::VectorX <prec>& solution,
                                         indexType numberBeamNodes);

  Types::Matrix33<prec> getWmatrix(const Types::Vector3 <prec>& a_vector);
  Types::Matrix33<prec> getMmatrix(const Types::Vector3 <prec>& r0mj_vector, 
                                   const Types::Vector3 <prec>& hmj_vector,
                                   const Types::Vector3<prec> &omega_rot);

  Types::MatrixXX <prec> getK_sigmaMatrix(const Types::VectorX <prec>& shapeDerivative,
                                          const Types::VectorX<prec>& shape,
                                          const Types::Matrix33<prec>& R01,
                                          const Types::Matrix33<prec>& R02,
                                          const Types::Matrix3X<prec>& coors,
                                          const Types::VectorX <prec>& solution,
                                          const Types::VectorX <prec>& sigma,
                                          indexType numDofs, indexType numberBeamNodes);


  void B_ei_ki_prime_matrix(const Types::Matrix33<prec>& r0I,
      const Types::Matrix33<prec>& r0,
      const Types::Matrix33<prec>& r0_prime,
      const Types::Vector3<prec>& xo_prime,
      Types::Matrix33<prec>* B_ei_mat,
      Types::Matrix33<prec>* B_ki_mat,
      Types::Matrix33<prec>* B_ki_prime_mat);

  void setTangentResidualPrismDispFormulation(
    PointerCollection& pointers,
    FiniteElement::LinearPrism *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs
  );

  void setTangentResidualPrismHuWashizuFormulation(
    PointerCollection& pointers,
    FiniteElement::LinearPrism* elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic>& stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1>& residual, std::vector<DegreeOfFreedom*>& Dofs);
  
  void setTangentResidualPrismgeononlinearDispFormulation(
    PointerCollection& pointers,
    FiniteElement::LinearPrism* elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic>& stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1>& residual, std::vector<DegreeOfFreedom*>& Dofs
  );
  void setTangentResidualPrismgeononlinearHuWashizuFormulation(
    PointerCollection& pointers,
    FiniteElement::LinearPrism* elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic>& stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1>& residual, std::vector<DegreeOfFreedom*>& Dofs
  );
  indexType meshIdDisp, meshIdRot, intOrderDisp, mode;
  prec EA,GA, GAy, GAz,GIo, EIy, EIz, EIyz, EIx;
  prec ys, zs;


  const static HistoryDataStructure m_historyDataStructure;
};
} // namespace HierAMuS
