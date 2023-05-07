// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "MatrixTypes.h"
#include "equations/DegreeOfFreedom.h"
#include <iostream>
#include <limits>


#include <elementFormulations/EL104_TimoshenkoPrism.h>
#include <elementFormulations/GenericElementFormulation.h>

#include <finiteElements/ElementTypes.h>
#include <finiteElements/GenericFiniteElement.h>
#include <finiteElements/LinearPrism.h>

#include <geometry/Base.h>
#include <geometry/Vertex.h>

#include <equations/EquationHandler.h>
#include <equations/GenericNodes.h>
#include <equations/NodeSet.h>

#include <solver/GenericSolutionState.h>

#include <math/geometry.h>

#include <pointercollection/pointercollection.h>

#include <Eigen/Dense>
#include <cmath>
#include <vector>

#include <limits>

#include "vtkCellType.h"

namespace HierAMuS::Elementformulations {

EL104_TimoshenkoPrism::EL104_TimoshenkoPrism(PointerCollection *ptrCol)
    : GenericElementFormulation(ptrCol) {}

EL104_TimoshenkoPrism::~EL104_TimoshenkoPrism() = default;

void EL104_TimoshenkoPrism::readData(PointerCollection &pointers,
                                     ParameterList &list) {
  std::string temp;
  this->meshIdDisp = list.getIndexVal("meshiddisp");
  this->meshIdRot = list.getIndexVal("meshidrot");
  this->intOrderDisp = list.getIndexVal("intorder");
  this->EA = list.getPrecVal("ea");
  this->GA = list.getPrecVal("ga");
  this->GIo = list.getPrecVal("gio");
  this->GAy = list.getPrecVal("gay");
  this->GAz = list.getPrecVal("gaz");
  this->EIy = list.getPrecVal("eiy");
  this->EIz = list.getPrecVal("eiz");
  this->EIx = list.getPrecVal("eix");
  this->EIyz = list.getPrecVal("eiyz");
  this->zs = list.getPrecVal("zs");
  this->ys = list.getPrecVal("ys");
  this->mode = list.getIndexVal("mode");
}

void EL104_TimoshenkoPrism::setDegreesOfFreedom(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) {
  switch (elem->getType()) {
  case FiniteElement::Elementtypes::LinearPrism: {
    FiniteElement::LinearPrism *el =
        static_cast<FiniteElement::LinearPrism *>(elem);
    el->setH1BeamShapes(pointers, this->meshIdDisp, this->intOrderDisp);
    el->setH1BeamShapes(pointers, this->meshIdRot, this->intOrderDisp);
    
    auto GP = el->getIntegrationPoints(pointers);
    
    GP.setTypeOrder(IntegrationType::Gauss1D, this->intOrderDisp);


  } break;
  default:
    throw std::runtime_error("Elementtype not supported");
  }
}

auto EL104_TimoshenkoPrism::getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> std::vector<DegreeOfFreedom *>
{

    FiniteElement::LinearPrism *passelem =
        static_cast<FiniteElement::LinearPrism *>(elem);
  std::vector<DegreeOfFreedom*> Dofs;
  std::vector<DegreeOfFreedom *> tempDofs;

  passelem->getH1BeamDofs(
      pointers, tempDofs, this->meshIdDisp,
      this->intOrderDisp); // gives degree of freedoms of displacemnets of nodes
                           // 1 and 2.[u1_x,u1_y,u1_z,u2_x,u2_y,u2_z]
  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());
  tempDofs.clear();
  passelem->getH1BeamDofs(
      pointers, tempDofs, this->meshIdRot,
      this->intOrderDisp); // gives degree of freedoms of rotations of nodes 1
                           // and 2.[w1_x,w1_y,w1_z,w2_x,w2_y,w2_z]
  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());
  return Dofs;
}

void EL104_TimoshenkoPrism::setTangentResidual(
  PointerCollection& pointers,
  FiniteElement::GenericFiniteElement *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  switch (elem->getType()) {
  case FiniteElement::Elementtypes::LinearPrism: {
    FiniteElement::LinearPrism *passelem =
        static_cast<FiniteElement::LinearPrism *>(elem);

    switch (this->mode) {
    case 1: // Displacement based formulation
      this->setTangentResidualPrismDispFormulation(pointers, passelem,
                                                     stiffness, residual, Dofs);
      break;
    case 2: // Hu-Washizu formulation
      this->setTangentResidualPrismHuWashizuFormulation(
          pointers, passelem,
                                                          stiffness, residual, Dofs);
      break;
    case 3: // Hu-Washizu Geometric non-linearity formulation
      this->setTangentResidualPrismgeononlinearDispFormulation(
          pointers, passelem, stiffness, residual, Dofs);
      break;
    case 4: // Hu-Washizu Geometric non-linearity formulation
      this->setTangentResidualPrismgeononlinearHuWashizuFormulation(
          pointers, passelem, stiffness, residual, Dofs);
      break;
    }
    break;
  }
  default:
    throw std::runtime_error(
        "Element type or it's mode is not compatible with element formulation");
  }
}

Types::Matrix66<prec> EL104_TimoshenkoPrism::getMaterialTangent() {

  Types::Matrix66<prec> matrix;
  matrix.setZero();

  matrix(0, 0) = this->EA;
  matrix(1, 1) = this->GAy;
  matrix(2, 2) = this->GAz;
  matrix(3, 3) = this->EIx;
  matrix(4, 4) = this->EIz;
  matrix(5, 5) = this->EIy;
  matrix(4, 5) = this->EIyz;
  matrix(5, 4) = this->EIyz;

  return matrix;
}

Types::Matrix66<prec> EL104_TimoshenkoPrism::getnonMaterialTangent() {

  Types::Matrix66<prec> matrix;
  matrix.setZero();

  matrix(0, 0) = this->EA;
  matrix(1, 1) = this->GA;
  matrix(2, 2) = this->GA;
  matrix(3, 3) = this->GIo;
  matrix(4, 4) = this->EIz;
  matrix(5, 5) = this->EIy;
  matrix(4, 5) = -(this->EIyz);
  matrix(5, 4) = -(this->EIyz);
  matrix(0, 4) = (this->EA) * (this->zs);
  matrix(0, 5) = -(this->EA) * (this->ys);
  matrix(4, 0) = (this->EA) * (this->zs);
  matrix(5, 0) = -(this->EA) * (this->ys);
  matrix(1, 3) = -(this->GAz);
  matrix(3, 1) = -(this->GAz);
  matrix(2, 3) = this->GAy;
  matrix(3, 2) = this->GAy;

  return matrix;
}

void EL104_TimoshenkoPrism::setTangentResidualPrismDispFormulation(
  PointerCollection& pointers,
  FiniteElement::LinearPrism *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  auto dmat = this->getMaterialTangent();

  auto GP = elem->getIntegrationPoints(pointers);
  GP.setTypeOrder(IntegrationType::Gauss1D, this->intOrderDisp);

  // R = [A1, A2, A3] as column vectors
  auto R01 = elem->getStartTriad(
      pointers); // orthonormal basis vector at node 1 of element
  auto R02 = elem->getEndTriad(
      pointers); // orthonormal basis vector at node 2 of element

  Types::VectorX<prec> shape, shapeDerivative;

  // Set up Degree of Freedom vector
  Dofs.clear();
  std::vector<DegreeOfFreedom *> tempDofs;

  elem->getH1BeamDofs(
      pointers, tempDofs, this->meshIdDisp,
      this->intOrderDisp); // gives degree of freedoms of displacemnets of nodes
                           // 1 and 2.[u1_x,u1_y,u1_z,u2_x,u2_y,u2_z]
  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());
  tempDofs.clear();
  elem->getH1BeamDofs(
      pointers, tempDofs, this->meshIdRot,
      this->intOrderDisp); // gives degree of freedoms of rotations of nodes 1
                           // and 2.[w1_x,w1_y,w1_z,w2_x,w2_y,w2_z]
  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());

  indexType numDofs = static_cast<indexType>(Dofs.size());

  // Set up size of stiffness matrix and residual vector.
  stiffness.resize(Dofs.size(), Dofs.size());
  residual.resize(Dofs.size());

  stiffness.setZero();
  residual.setZero();

  Types::MatrixXX<prec> BMatrix;
  BMatrix.resize(6, Dofs.size());
  BMatrix.setZero();

  indexType numberBeamNodes = Dofs.size() / indexType(6);

  // Get the solution vector
  Types::VectorX<prec> solution;
  solution = pointers.getSolutionState()->getSolution(
      Dofs); // Displacements and rotations at x,y,z directions of nodes of
             // element

  for (auto i = 0; i < GP.getTotalGP(); ++i) {
    prec jacobian = elem->getBeamJacobian(pointers, GP.getXi(i));
    elem->getH1BeamShapes(pointers, this->intOrderDisp, jacobian, shape,
                          shapeDerivative, GP.getXi(i));

    auto R0 =
        shape(0) * R01 + shape(1) * R02; // orthonormal basis vector of element

    BMatrix =
        this->getBMatrix(shapeDerivative, shape, R0, numDofs, numberBeamNodes);

    stiffness +=
        BMatrix.transpose() * dmat * BMatrix * jacobian * GP.getWeight(i);
  }
  residual = stiffness * solution;

  //  Eigen::EigenSolver<Types::MatrixXX<prec>> es(stiffness);
  //  std::cout << "Eigenvalues of J are: " << es.eigenvalues() << std::endl;
}

Types::Matrix6X<prec> EL104_TimoshenkoPrism::getBMatrix(
  const Types::VectorX<prec> &shapeDerivative,
  const Types::VectorX<prec> &shape, const Types::Matrix33<prec> &R0,
  indexType numDofs, indexType numberBeamNodes) {

  Types::Vector3<prec> A01;
  Types::Vector3<prec> A02;
  Types::Vector3<prec> A03;
  Types::Matrix6X<prec> BMatrix;

  A01 = R0.block(0, 0, 3, 1);
  A02 = R0.block(0, 1, 3, 1);
  A03 = R0.block(0, 2, 3, 1);

  BMatrix.resize(6, numDofs);
  BMatrix.setZero();

  for (auto j = 0; j < numberBeamNodes; ++j) {

    BMatrix.block(0, 3 * j, 1, 3) =
        A01.transpose() * shapeDerivative(j); // (epsilon x), epsilon xi1
    BMatrix.block(1, 3 * j, 1, 3) =
        A02.transpose() *
        shapeDerivative(j); // gamma y, gamma xi2; from displacements
    BMatrix.block(2, 3 * j, 1, 3) =
        A03.transpose() *
        shapeDerivative(j); // gamma z, gamma xi3; from displacements

    indexType offset = 3 * numberBeamNodes;
    BMatrix.block(1, 3 * j + offset, 1, 3) =
        -A03.transpose() * shape(j); // gamma y, gamma xi2; from rotations
    BMatrix.block(2, 3 * j + offset, 1, 3) =
        A02.transpose() * shape(j); // gamma z, gamma xi3; from rotations

    BMatrix.block(3, 3 * j + offset, 1, 3) =
        A01.transpose() * shapeDerivative(j); // kappa x, kappa xi1
    BMatrix.block(4, 3 * j + offset, 1, 3) =
        A02.transpose() * shapeDerivative(j);
    BMatrix.block(5, 3 * j + offset, 1, 3) =
        -A03.transpose() * shapeDerivative(j);
  }
  return BMatrix;
}

void EL104_TimoshenkoPrism::setTangentResidualPrismHuWashizuFormulation(
  PointerCollection& pointers,
  FiniteElement::LinearPrism *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  auto dmat = this->getMaterialTangent();

  auto GP = pointers.getIntegrationPoints(-1);
  GP.setTypeOrder(IntegrationType::Gauss1D, this->intOrderDisp);

  // R = [A1, A2, A3] as column vectors
  auto R01 = elem->getStartTriad(
      pointers); // orthonormal basis vector at node 1 of element
  auto R02 = elem->getEndTriad(
      pointers); // orthonormal basis vector at node 2 of element

  Types::VectorX<prec> shape, shapeDerivative;

  // Set up Degree of Freedom vector
  Dofs.clear();
  std::vector<DegreeOfFreedom *> tempDofs;

  elem->getH1BeamDofs(
      pointers, tempDofs, this->meshIdDisp,
      this->intOrderDisp); // gives degree of freedoms of displacemnets of nodes
                           // 1 and 2.[u1_x,u1_y,u1_z,u2_x,u2_y,u2_z]
  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());
  tempDofs.clear();
  elem->getH1BeamDofs(
      pointers, tempDofs, this->meshIdRot,
      this->intOrderDisp); // gives degree of freedoms of rotations of nodes 1
                           // and 2.[w1_x,w1_y,w1_z,w2_x,w2_y,w2_z]

  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());

  indexType numDofs = static_cast<indexType>(Dofs.size());

  // Set up size of stiffness matrix and residual vector.
  stiffness.resize(Dofs.size(), Dofs.size());
  residual.resize(Dofs.size());

  stiffness.setZero();
  residual.setZero();

  Types::MatrixXX<prec> BMatrix;
  BMatrix.resize(6, numDofs);
  BMatrix.setZero();

  Types::Matrix6X<prec> Amat;
  Amat.resize(6, 6);

  Amat = this->getLocalStrainStressInterpolation();

  indexType numberBeamNodes = Dofs.size() / indexType(6);

  Types::MatrixXX<prec> G, H, L;
  G.resize(Dofs.size(), 6);
  H.resize(6, 6);
  L.resize(6, 6);

  G.setZero();
  H.setZero();
  L.setZero();

  // Get the solution vector
  Types::VectorX<prec> solution;
  solution = pointers.getSolutionState()->getSolution(
      Dofs); // Displacements and rotations at x,y,z directions of nodes of
             // element

  for (auto i = 0; i < GP.getTotalGP(); ++i) {
    prec jacobian = elem->getBeamJacobian(pointers, GP.getXi(i));
    elem->getH1BeamShapes(pointers, this->intOrderDisp, jacobian, shape,
                          shapeDerivative, GP.getXi(i));

    auto R0 =
        shape(0) * R01 + shape(1) * R02; // orthonormal basis vector of element

    BMatrix =
        this->getBMatrix(shapeDerivative, shape, R0, numDofs, numberBeamNodes);

    G += BMatrix.transpose() * Amat * jacobian * GP.getWeight(i);
    H -= Amat.transpose() * Amat * jacobian * GP.getWeight(i);
    L += Amat.transpose() * dmat * Amat * jacobian * GP.getWeight(i);
  }

  stiffness = G * (H.transpose()).inverse() * L * H.inverse() * G.transpose();
  residual = stiffness * solution;
}

Types::Matrix3X<prec>
EL104_TimoshenkoPrism::Hmatrix(const Types::VectorX<prec> &solution,
                               indexType numberBeamNodes) {

  Types::Matrix33<prec> omega, Identity, Hmat_I;
  Types::Vector3<prec> omega_rot; // vector of the rotations at each node

  Types::Matrix3X<prec> Hmat;
  Hmat.resize(3, 3 * numberBeamNodes);
  Hmat.setZero();
  Hmat_I.setZero();
  omega.setZero();

  Identity.resize(3, 3);
  Identity.setZero();

  for (auto i = 0; i < 3; ++i) {
    Identity(i, i) = prec(1);
  }

  for (auto j = 0; j < numberBeamNodes; ++j) {

    omega_rot = solution.segment(
        3 * j + 6,
        3); // comes out from the solution matrix after the displacements
    prec omega_rot_norm =
        omega_rot.norm(); // gives norm of the rotation of the nodes

    if (omega_rot_norm <= std::numeric_limits<prec>::epsilon() * 100) {
      Hmat_I = Identity;
    } else {
      omega(0, 1) = -omega_rot[2];
      omega(0, 2) = omega_rot[1];
      omega(1, 0) = omega_rot[2];
      omega(1, 2) = -omega_rot[0];
      omega(2, 0) = -omega_rot[1];
      omega(2, 1) = omega_rot[0];

      Hmat_I =
          Identity +
          (((prec)1 - cos(omega_rot_norm)) / omega_rot_norm / omega_rot_norm) *
              omega +
          ((omega_rot_norm - sin(omega_rot_norm)) / omega_rot_norm /
           omega_rot_norm / omega_rot_norm) *
              omega * omega;
    }
    Hmat.block(0, 3 * j, 3, 3) = Hmat_I;
  }

  return Hmat;
}

Types::Matrix3X<prec>
EL104_TimoshenkoPrism::Rmatrix(const Types::VectorX<prec> &solution,
                               indexType numberBeamNodes) {

  /* Types::Matrix3X<prec> Rmat;*/
  Types::Matrix33<prec> omega, Identity, Rmat_I;
  Types::Vector3<prec> omega_rot; // vector of the rotations at each node

  Types::Matrix3X<prec> Rmat;
  Rmat.resize(3, 3 * numberBeamNodes);
  Rmat.setZero();
  Rmat_I.setZero();
  omega.setZero();

  Identity.resize(3, 3);
  Identity.setZero();

  for (auto i = 0; i < 3; ++i) {
    Identity(i, i) = prec(1);
  }

  for (auto j = 0; j < numberBeamNodes; ++j) {

    omega_rot = solution.segment(
        3 * j + 6,
        3); // comes out from the solution matrix after the displacements
    prec omega_rot_norm =
        omega_rot.norm(); // gives norm of the rotation of the nodes
    if (omega_rot_norm <= std::numeric_limits<prec>::epsilon() * 100) {
      Rmat_I = Identity;
    } else {
      omega(0, 1) = -omega_rot[2];
      omega(0, 2) = omega_rot[1];
      omega(1, 0) = omega_rot[2];
      omega(1, 2) = -omega_rot[0];
      omega(2, 0) = -omega_rot[1];
      omega(2, 1) = omega_rot[0];

      Rmat_I =
          Identity + (sin(omega_rot_norm) / omega_rot_norm) * omega +
          (((prec)1 - cos(omega_rot_norm)) / omega_rot_norm / omega_rot_norm) *
              omega * omega;
    }
    Rmat.block(0, 3 * j, 3, 3) = Rmat_I;
  }
  return Rmat;
}

Types::Vector3<prec> EL104_TimoshenkoPrism::xo_prime_matrix(
  const Types::VectorX<prec> &shapeDerivative,
  const Types::VectorX<prec> &solution, const Types::Matrix3X<prec>& coors,
  indexType numberBeamNodes) {

  Types::Vector3<prec> xo_prime, u_disp;
  xo_prime.setZero();

  for (auto j = 0; j < numberBeamNodes; ++j) {

    u_disp = solution.segment(3 * j, 3);

    xo_prime += shapeDerivative(j) * (coors.block(0, j, 3, 1) + u_disp);
  }

  return xo_prime;
}

Types::Vector3<prec> EL104_TimoshenkoPrism::Xo_prime_matrix(
  const Types::VectorX<prec> &shapeDerivative,
  const Types::Matrix3X<prec>& coors, indexType numberBeamNodes) {

  Types::Vector3<prec> Xo_prime;
  Xo_prime.setZero();

  for (auto j = 0; j < numberBeamNodes; ++j) {

    Xo_prime += shapeDerivative(j) * coors.block(0, j, 3, 1);
  }

  return Xo_prime;
}

void EL104_TimoshenkoPrism::B_ei_ki_prime_matrix(
    const Types::Matrix33<prec> &r0I, const Types::Matrix33<prec> &r0,
    const Types::Matrix33<prec> &r0_prime, const Types::Vector3<prec> &xo_prime,
    Types::Matrix33<prec> *B_ei_mat, Types::Matrix33<prec> *B_ki_mat,
    Types::Matrix33<prec> *B_ki_prime_mat) {

  Types::Matrix33<prec> B_ei, B_ki, B_ki_prime;

  Types::Vector3<prec> a01, a02, a03, aI1, aI2, aI3, a01_prime, a02_prime,
      a03_prime;

  a01 = r0.block(0, 0, 3, 1);
  a02 = r0.block(0, 1, 3, 1);
  a03 = r0.block(0, 2, 3, 1);

  a01_prime = r0_prime.block(0, 0, 3, 1);
  a02_prime = r0_prime.block(0, 1, 3, 1);
  a03_prime = r0_prime.block(0, 2, 3, 1);

  aI1 = r0I.block(0, 0, 3, 1);
  aI2 = r0I.block(0, 1, 3, 1);
  aI3 = r0I.block(0, 2, 3, 1);

  B_ei.block(0, 0, 3, 1) = aI1.cross(xo_prime);
  B_ei.block(0, 1, 3, 1) = aI2.cross(xo_prime);
  B_ei.block(0, 2, 3, 1) = aI3.cross(xo_prime);

  B_ki.block(0, 0, 3, 1) = aI2.cross(a03);
  B_ki.block(0, 1, 3, 1) = aI3.cross(a01);
  B_ki.block(0, 2, 3, 1) = aI1.cross(a02);

  B_ki_prime.block(0, 0, 3, 1) = aI3.cross(a02_prime);
  B_ki_prime.block(0, 1, 3, 1) = aI1.cross(a03_prime);
  B_ki_prime.block(0, 2, 3, 1) = aI2.cross(a01_prime);

  *B_ei_mat = B_ei;
  *B_ki_mat = B_ki;
  *B_ki_prime_mat = B_ki_prime;
};

Types::Matrix6X<prec> EL104_TimoshenkoPrism::getNonBMatrix(
  const Types::VectorX<prec> &shapeDerivative,
  const Types::VectorX<prec> &shape, const Types::Matrix33<prec> &R01,
  const Types::Matrix33<prec> &R02, const Types::Matrix3X<prec> &coors,
  const Types::VectorX<prec> &solution, indexType numDofs,
  indexType numberBeamNodes) {

  Types::Vector3<prec> A01, A02, A03;
  Types::Matrix6X<prec> BMatrix;

  auto R0 = shape(0) * R01 + shape(1) * R02;

  Types::Matrix3X<prec> R0I;
  R0I.resize(3, 3 * numberBeamNodes);
  R0I.block(0, 0, 3, 3) = R01;
  R0I.block(0, 3, 3, 3) = R02;

  A01 = R0.block(0, 0, 3, 1);
  A02 = R0.block(0, 1, 3, 1);
  A03 = R0.block(0, 2, 3, 1);

  BMatrix.resize(6, numDofs);
  BMatrix.setZero();

  Types::Matrix3X<prec> Rmat, Hmat;
  Rmat.resize(3, 3 * numberBeamNodes);
  Rmat.setZero();

  Hmat.resize(3, 3 * numberBeamNodes);
  Hmat.setZero();

  Rmat = this->Rmatrix(solution, numberBeamNodes);

  Hmat = this->Hmatrix(solution, numberBeamNodes);

  Types::Matrix33<prec> r01, r02, r0, r0_prime, r0I;
  r01 = Rmat.block(0, 0, 3, 3) * R01;
  r02 = Rmat.block(0, 3, 3, 3) * R02;

  Types::Matrix3X<prec> RR;
  RR.resize(3, 6);
  RR.block(0, 0, 3, 3) = Rmat.block(0, 0, 3, 3) * R01;
  RR.block(0, 3, 3, 3) = Rmat.block(0, 3, 3, 3) * R02;

  r0 = shape(0) * r01 + shape(1) * r02;

  r0_prime = shapeDerivative(0) * r01 + shapeDerivative(1) * r02;

  Types::Vector3<prec> xo_prime;
  xo_prime =
      this->xo_prime_matrix(shapeDerivative, solution, coors, numberBeamNodes);

  for (auto j = 0; j < numberBeamNodes; ++j) {

    Types::Matrix33<prec> B_ei_mat, B_ki_mat, B_ki_prime_mat;
    B_ei_mat.setZero();
    B_ki_mat.setZero();
    B_ki_prime_mat.setZero();

    r0I = Rmat.block(0, 3 * j, 3, 3) * R0I.block(0, 3 * j, 3, 3);

    this->B_ei_ki_prime_matrix(r0I, r0, r0_prime, xo_prime, &B_ei_mat,
                               &B_ki_mat, &B_ki_prime_mat);

    BMatrix.block(0, 3 * j, 3, 3) = shapeDerivative(j) * r0.transpose();

    indexType offset = 3 * numberBeamNodes;

    BMatrix.block(0, 3 * j + offset, 3, 3) =
        shape(j) * B_ei_mat.transpose() * Hmat.block(0, 3 * j, 3, 3);

    BMatrix.block(3, 3 * j + offset, 3, 3) =
        (shapeDerivative(j) * B_ki_mat.transpose() +
         shape(j) * B_ki_prime_mat.transpose()) *
        Hmat.block(0, 3 * j, 3, 3);
  }

  return BMatrix;
};

Types::Vector6<prec> EL104_TimoshenkoPrism::getstrain_Vector(
  const Types::VectorX<prec> &shapeDerivative,
  const Types::VectorX<prec> &shape, const Types::Matrix33<prec> &R01,
  const Types::Matrix33<prec> &R02, const Types::Matrix3X<prec> &coors,
  const Types::VectorX<prec> &solution, indexType numberBeamNodes) {

  Types::Vector6<prec> strain_vector;
  Types::Vector3<prec> A01, A02, A03, a01, a02, a03, A01_prime, A02_prime,
      A03_prime, a01_prime, a02_prime, a03_prime;

  auto R0 = shape(0) * R01 + shape(1) * R02;
  auto R0_prime = shapeDerivative(0) * R01 + shapeDerivative(1) * R02;

  A01 = R0.block(0, 0, 3, 1);
  A02 = R0.block(0, 1, 3, 1);
  A03 = R0.block(0, 2, 3, 1);

  A01_prime = R0_prime.block(0, 0, 3, 1);
  A02_prime = R0_prime.block(0, 1, 3, 1);
  A03_prime = R0_prime.block(0, 2, 3, 1);

  Types::Matrix3X<prec> Rmat;
  Rmat.resize(3, 3 * numberBeamNodes);
  Rmat.setZero();

  Rmat = this->Rmatrix(solution, numberBeamNodes);

  Types::Matrix33<prec> r01, r02, r0, r0_prime;
  r01 = Rmat.block(0, 0, 3, 3) * R01;
  r02 = Rmat.block(0, 3, 3, 3) * R02;
  // std::cout << R01 << std::endl;
  // std::cout << R02 << std::endl;

  r0 = shape(0) * r01 + shape(1) * r02;

  r0_prime = shapeDerivative(0) * r01 + shapeDerivative(1) * r02;

  a01 = r0.block(0, 0, 3, 1);
  a02 = r0.block(0, 1, 3, 1);
  a03 = r0.block(0, 2, 3, 1);

  a01_prime = r0_prime.block(0, 0, 3, 1);
  a02_prime = r0_prime.block(0, 1, 3, 1);
  a03_prime = r0_prime.block(0, 2, 3, 1);

  Types::Vector3<prec> xo_prime, Xo_prime;
  xo_prime =
      this->xo_prime_matrix(shapeDerivative, solution, coors, numberBeamNodes);
  Xo_prime = this->Xo_prime_matrix(shapeDerivative, coors, numberBeamNodes);

  strain_vector(0) = a01.dot(xo_prime) - A01.dot(Xo_prime);
  strain_vector(1) = a02.dot(xo_prime) - A02.dot(Xo_prime);
  strain_vector(2) = a03.dot(xo_prime) - A03.dot(Xo_prime);
  strain_vector(3) = a03.dot(a02_prime) - A03.dot(A02_prime);
  strain_vector(4) = a01.dot(a03_prime) - A01.dot(A03_prime);
  strain_vector(5) = a02.dot(a01_prime) - A02.dot(A01_prime);

  return strain_vector;
};

Types::Matrix33<prec>
EL104_TimoshenkoPrism::getWmatrix(const Types::Vector3<prec> &a_vector) {

  Types::Matrix33<prec> Wmatrix;
  Wmatrix.setZero();

  Wmatrix(0, 1) = -a_vector(2);
  Wmatrix(0, 2) = a_vector(1);
  Wmatrix(1, 0) = a_vector(2);
  Wmatrix(1, 2) = -a_vector(0);
  Wmatrix(2, 0) = -a_vector(1);
  Wmatrix(2, 1) = a_vector(0);

  return Wmatrix;
};

Types::Matrix33<prec>
EL104_TimoshenkoPrism::getMmatrix(const Types::Vector3<prec> &r0mj_vector,
                                  const Types::Vector3<prec> &hmj_vector,
                                  const Types::Vector3<prec> &omega_rot) {

  Types::Matrix33<prec> M_mjmatrix;
  M_mjmatrix.setZero();

  prec omega_rot_norm = omega_rot.norm();

  if (omega_rot_norm <= std::numeric_limits<prec>::epsilon() * 100) {

    M_mjmatrix = (prec)0.5 * r0mj_vector * hmj_vector.transpose() +
                 (prec)0.5 * hmj_vector * r0mj_vector.transpose();
  } else {
    Types::Matrix33<prec> Identity;

    Identity.setZero();

    for (auto i = 0; i < 3; ++i) {
      Identity(i, i) = prec(1);
    }

    Types::Vector3<prec> b_mj, r0mj_hat;
    b_mj.setZero();
    r0mj_hat.setZero();

    b_mj = r0mj_vector.cross(hmj_vector);

    prec temp = b_mj.dot(omega_rot);
    prec temp_1 = r0mj_vector.dot(hmj_vector);

    // prec c_1 = ((prec)1 - cos(omega_rot_norm)) / omega_rot_norm /
    // omega_rot_norm;
    prec c_2 = (omega_rot_norm - sin(omega_rot_norm)) / omega_rot_norm /
               omega_rot_norm / omega_rot_norm;
    prec tempVal = cos(omega_rot_norm) - prec(1);
    prec c_3;
    if (tempVal <= std::numeric_limits<prec>::epsilon() * 100) {
      c_3 = prec(1) / prec(6);
    } else {
      c_3 = (omega_rot_norm * sin(omega_rot_norm) +
             (prec)2 * cos(omega_rot_norm) - (prec)2) /
            (omega_rot_norm * omega_rot_norm * (cos(omega_rot_norm) - (prec)1));
    }

    // prec c_4 = -c_1 * c_3;
    prec c_6 = (c_3 - c_2) / omega_rot_norm / omega_rot_norm;
    // prec c_5 = c_6 - c_2 * c_3;
    // prec c_7 = (prec)1 - (prec)0.5 * c_3* omega_rot_norm * omega_rot_norm;
    // prec c_8 = (prec)0.25 - (prec)0.5 * c_3 + (prec)0.25 * c_3 * c_3 *
    // omega_rot_norm * omega_rot_norm;
    prec c_9 = -(prec)0.25 + c_3 -
               (prec)0.25 * c_3 * c_3 * omega_rot_norm * omega_rot_norm;
    prec c_10 =
        c_2 * ((prec)1 - c_9 * omega_rot_norm * omega_rot_norm) * temp - temp_1;

    r0mj_hat = -c_3 * b_mj + (c_6 + c_2 * c_9) * temp * omega_rot;

    M_mjmatrix = (prec)0.5 * r0mj_vector * hmj_vector.transpose() +
                 (prec)0.5 * hmj_vector * r0mj_vector.transpose() +
                 (prec)0.5 * r0mj_hat * omega_rot.transpose() +
                 (prec)0.5 * omega_rot * r0mj_hat.transpose() + c_10 * Identity;
  }

  return M_mjmatrix;
};

Types::MatrixXX<prec> EL104_TimoshenkoPrism::getK_sigmaMatrix(
  const Types::VectorX<prec> &shapeDerivative,
  const Types::VectorX<prec> &shape, const Types::Matrix33<prec> &R01,
  const Types::Matrix33<prec> &R02, const Types::Matrix3X<prec> &coors,
  const Types::VectorX<prec> &solution, const Types::VectorX<prec> &sigma,
  indexType numDofs, indexType numberBeamNodes) {

  Types::Vector3<prec> A01, A02, A03;

  auto R0 = shape(0) * R01 + shape(1) * R02;

  Types::Matrix3X<prec> R0I;
  R0I.resize(3, 3 * numberBeamNodes);
  R0I.block(0, 0, 3, 3) = R01;
  R0I.block(0, 3, 3, 3) = R02;

  A01 = R0.block(0, 0, 3, 1);
  A02 = R0.block(0, 1, 3, 1);
  A03 = R0.block(0, 2, 3, 1);

  Types::MatrixXX<prec> K_sigmaMatrix;
  K_sigmaMatrix.resize(numDofs, numDofs);
  K_sigmaMatrix.setZero();

  Types::Matrix3X<prec> Rmat, Hmat;
  Rmat.resize(3, 3 * numberBeamNodes);
  Rmat.setZero();

  Hmat.resize(3, 3 * numberBeamNodes);
  Hmat.setZero();

  Rmat = this->Rmatrix(solution, numberBeamNodes);

  Hmat = this->Hmatrix(solution, numberBeamNodes);

  Types::Matrix33<prec> r01, r02, r0, r0_prime, r0j, r0k;
  r01 = Rmat.block(0, 0, 3, 3) * R01;
  r02 = Rmat.block(0, 3, 3, 3) * R02;

  r0 = shape(0) * r01 + shape(1) * r02;

  r0_prime = shapeDerivative(0) * r01 + shapeDerivative(1) * r02;

  Types::Vector3<prec> x0_prime;
  x0_prime =
      this->xo_prime_matrix(shapeDerivative, solution, coors, numberBeamNodes);

  for (auto j = 0; j < numberBeamNodes; ++j) {

    for (auto k = 0; k < numberBeamNodes; ++k) {

      Types::Matrix33<prec> Wjk_1, Wjk_2;
      Wjk_1.setZero();
      Wjk_2.setZero();

      r0j = Rmat.block(0, 3 * j, 3, 3) * R0I.block(0, 3 * j, 3, 3);
      r0k = Rmat.block(0, 3 * k, 3, 3) * R0I.block(0, 3 * k, 3, 3);

      auto r01j = r0j.block(0, 0, 3, 1);
      auto r02j = r0j.block(0, 1, 3, 1);
      auto r03j = r0j.block(0, 2, 3, 1);

      auto W1j = this->getWmatrix(r01j);
      auto W2j = this->getWmatrix(r02j);
      auto W3j = this->getWmatrix(r03j);

      auto r01k = r0k.block(0, 0, 3, 1);
      auto r02k = r0k.block(0, 1, 3, 1);
      auto r03k = r0k.block(0, 2, 3, 1);

      auto W1k = this->getWmatrix(r01k);
      auto W2k = this->getWmatrix(r02k);
      auto W3k = this->getWmatrix(r03k);

      Wjk_1 = sigma(3) * W2j * W3k.transpose() +
              sigma(4) * W3j * W1k.transpose() +
              sigma(5) * W1j * W2k.transpose();
      Wjk_2 = sigma(3) * W3j * W2k.transpose() +
              sigma(4) * W1j * W3k.transpose() +
              sigma(5) * W2j * W1k.transpose();

      Types::Vector3<prec> F_r0j, F_r0k;
      F_r0j = sigma(0) * r01j + sigma(1) * r02j + sigma(2) * r03j;
      auto Wfj = this->getWmatrix(F_r0j);

      F_r0k = sigma(0) * r01k + sigma(1) * r02k + sigma(2) * r03k;
      auto Wfk = this->getWmatrix(F_r0k);

      Types::Vector3<prec> h_1j, h_2j, h_3j;

      h_1j = shape(j) * sigma(0) * x0_prime +
             shape(j) * sigma(4) * r0_prime.block(0, 2, 3, 1) +
             shapeDerivative(j) * sigma(5) * r0.block(0, 1, 3, 1);
      h_2j = shape(j) * sigma(1) * x0_prime +
             shape(j) * sigma(5) * r0_prime.block(0, 0, 3, 1) +
             shapeDerivative(j) * sigma(3) * r0.block(0, 2, 3, 1);
      h_3j = shape(j) * sigma(2) * x0_prime +
             shape(j) * sigma(3) * r0_prime.block(0, 1, 3, 1) +
             shapeDerivative(j) * sigma(4) * r0.block(0, 0, 3, 1);

      Types::Matrix33<prec> MMatrix, M_1J, M_2J, M_3J;
      MMatrix.setZero();
      M_1J.setZero();
      M_2J.setZero();
      M_3J.setZero();

      if (j == k) {
        // Types::Matrix33<prec>MMatrix_1I, MMatrix_2I, MMatrix_3I;
        Types::Vector3<prec> omega_rot;
        omega_rot = solution.segment(3 * j + 6, 3);

        M_1J = this->getMmatrix(r01j, h_1j, omega_rot);
        M_2J = this->getMmatrix(r02j, h_2j, omega_rot);
        M_3J = this->getMmatrix(r03j, h_3j, omega_rot);

        // MMatrix = this->getMmatrix(r01j, h_1j,omega_rot) +
        // this->getMmatrix(r02j, h_2j, omega_rot) + this->getMmatrix(r03j,
        // h_3j, omega_rot);

        MMatrix = M_1J + M_2J + M_3J;
      }
      // std::cout << MMatrix << std::endl;
      Types::Matrix33<prec> G_jkMatrix;
      G_jkMatrix.setZero();

      G_jkMatrix = shapeDerivative(j) * shape(k) * Wjk_1 +
                   shape(j) * shapeDerivative(k) * Wjk_2 + MMatrix;

      auto Hj = Hmat.block(0, 3 * j, 3, 3);
      auto Hk = Hmat.block(0, 3 * k, 3, 3);

      indexType offset = 3 * numberBeamNodes;

      K_sigmaMatrix.block(3 * j, 3 * k + offset, 3, 3) =
          shapeDerivative(j) * shape(k) * Wfk.transpose() * Hk;

      K_sigmaMatrix.block(3 * j + offset, 3 * k, 3, 3) =
          shape(j) * shapeDerivative(k) * Hj.transpose() * Wfj;

      K_sigmaMatrix.block(3 * j + offset, 3 * k + offset, 3, 3) =
          Hj.transpose() * G_jkMatrix * Hk;
    }
  }

  return K_sigmaMatrix;
};

void EL104_TimoshenkoPrism::setTangentResidualPrismgeononlinearDispFormulation(
  PointerCollection& pointers,
  FiniteElement::LinearPrism *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  auto dmat = this->getMaterialTangent();

  auto GP = pointers.getIntegrationPoints(-1);
  GP.setTypeOrder(IntegrationType::Gauss1D, this->intOrderDisp);

  // R = [A1, A2, A3] as column vectors
  auto R01 = elem->getStartTriad(pointers);
  auto R02 = elem->getEndTriad(pointers);

  Types::VectorX<prec> shape, shapeDerivative;

  // Set up Degree of Freedom vector
  Dofs.clear();
  std::vector<DegreeOfFreedom *> tempDofs;

  elem->getH1BeamDofs(pointers, tempDofs, this->meshIdDisp,
                      this->intOrderDisp);
  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());
  tempDofs.clear();
  elem->getH1BeamDofs(pointers, tempDofs, this->meshIdRot,
                      this->intOrderDisp);
  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());

  indexType numDofs = static_cast<indexType>(Dofs.size());

  // Set up size of stiffness matrix and residual vector.
  stiffness.resize(Dofs.size(), Dofs.size());
  residual.resize(Dofs.size());

  stiffness.setZero();
  residual.setZero();

  Types::MatrixXX<prec> BMatrix, K_sigmaMatrix;
  BMatrix.resize(6, Dofs.size());
  BMatrix.setZero();
  K_sigmaMatrix.resize(Dofs.size(), Dofs.size());
  K_sigmaMatrix.setZero();

  indexType numberBeamNodes = Dofs.size() / indexType(6);

  // Get the solution vector
  Types::VectorX<prec> solution, newtonSol;
  solution = pointers.getSolutionState()->getSolution(Dofs);
  newtonSol = pointers.getSolutionState()->getNewtonSolution(Dofs);

  // Strain Vector defining
  Types::Vector6<prec> strains, sigma;

  for (auto i = 0; i < GP.getTotalGP(); ++i) {
    prec jacobian = elem->getBeamJacobian(pointers, GP.getXi(i));
    elem->getH1BeamShapes(pointers, this->intOrderDisp, jacobian, shape,
                          shapeDerivative, GP.getXi(i));

    auto &vert1 = elem->getVertex(pointers, 0);
    auto &vert2 = elem->getVertex(pointers, 3);
    Types::Matrix3X<prec> coors;
    coors.resize(3, 2);
    coors.block(0, 0, 3, 1) = vert1.getCoordinates();
    coors.block(0, 1, 3, 1) = vert2.getCoordinates();

    BMatrix = this->getNonBMatrix(shapeDerivative, shape, R01, R02, coors,
                                  solution, numDofs, numberBeamNodes);

    strains = this->getstrain_Vector(shapeDerivative, shape, R01, R02, coors,
                                     solution, numberBeamNodes);
    sigma = dmat * strains;

    K_sigmaMatrix =
        this->getK_sigmaMatrix(shapeDerivative, shape, R01, R02, coors,
                               solution, sigma, numDofs, numberBeamNodes);

    stiffness +=
        BMatrix.transpose() * dmat * BMatrix * jacobian * GP.getWeight(i) +
        K_sigmaMatrix * jacobian * GP.getWeight(i);

    residual += BMatrix.transpose() * sigma * jacobian * GP.getWeight(i);
  }
};

void EL104_TimoshenkoPrism::
    setTangentResidualPrismgeononlinearHuWashizuFormulation(
  PointerCollection& pointers,
  FiniteElement::LinearPrism *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {
  auto dmat = this->getMaterialTangent();

  auto GP = pointers.getIntegrationPoints(-1);
  GP.setTypeOrder(IntegrationType::Gauss1D, this->intOrderDisp);


  auto histIt = elem->getHistoryDataIterator(pointers);
  auto sol = pointers.getSolutionState();
  auto G_previous = histIt.getFieldElementConst(0);
  auto F_previous = histIt.getFieldElementConst(1);
  auto H_previous = histIt.getFieldElementConst(2);
  auto f_s_previous = histIt.getFieldElementConst(3);
  auto f_e_previous = histIt.getFieldElementConst(4);
  auto epsilon_hat_previous = histIt.getFieldElementConst(5);
  auto sigma_hat_previous = histIt.getFieldElementConst(6);


  // R = [A1, A2, A3] as column vectors
  auto R01 = elem->getStartTriad(pointers);
  auto R02 = elem->getEndTriad(pointers);

  Types::VectorX<prec> shape, shapeDerivative;

  // Set up Degree of Freedom vector
  Dofs.clear();
  std::vector<DegreeOfFreedom *> tempDofs;

  elem->getH1BeamDofs(pointers, tempDofs, this->meshIdDisp,
                      this->intOrderDisp);
  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());
  tempDofs.clear();
  elem->getH1BeamDofs(pointers, tempDofs, this->meshIdRot,
                      this->intOrderDisp);
  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());

  indexType numDofs = static_cast<indexType>(Dofs.size());

  // Set up size of stiffness matrix and residual vector.
  stiffness.resize(Dofs.size(), Dofs.size());
  residual.resize(Dofs.size());

  stiffness.setZero();
  residual.setZero();

  Types::MatrixXX<prec> BMatrix, K_sigmaMatrix;
  BMatrix.resize(6, Dofs.size());
  BMatrix.setZero();
  K_sigmaMatrix.resize(Dofs.size(), Dofs.size());
  K_sigmaMatrix.setZero();

  indexType numberBeamNodes = Dofs.size() / indexType(6);

  // Get the solution vector
  Types::VectorX<prec> solution, newtonSol;
  solution = pointers.getSolutionState()->getSolution(Dofs);
  newtonSol = pointers.getSolutionState()->getNewtonSolution(Dofs);

  Types::Vector6<prec> epsilon_hat, sigma_hat, d_epsilon_hat, d_sigma_hat;
  epsilon_hat.setZero();
  sigma_hat.setZero();
  d_epsilon_hat.setZero();
  d_sigma_hat.setZero();

  Types::Matrix66<prec> nullMatrix;
  nullMatrix.setZero();

  if (F_previous != nullMatrix) {
    d_epsilon_hat = -(F_previous.inverse()).transpose() *
                    (G_previous * newtonSol + f_s_previous);
    d_sigma_hat =
        -F_previous.inverse() * (H_previous * d_epsilon_hat + f_e_previous);
  }

  epsilon_hat_previous += d_epsilon_hat;
  sigma_hat_previous += d_sigma_hat;

  epsilon_hat = epsilon_hat_previous;
  sigma_hat = sigma_hat_previous;

  // std::cout << newtonSol << std::endl;

  // Strain Vector defining
  Types::Vector6<prec> strains, sigma, strains_g;

  Types::Matrix66<prec> Amat;

  Amat = this->getLocalStrainStressInterpolation();

  Types::MatrixXX<prec> Kg, G, H, F;
  Kg.resize(Dofs.size(), Dofs.size());
  // G.resize(Dofs.size(),6);
  G.resize(6, Dofs.size());
  H.resize(6, 6);
  F.resize(6, 6);

  Kg.setZero();
  G.setZero();
  H.setZero();
  F.setZero();

  Types::VectorX<prec> f_i;
  f_i.resize(12);
  f_i.setZero();

  Types::Vector6<prec> f_e, f_e_e, f_s, f_s_e;
  f_e.setZero();
  f_e_e.setZero();
  f_s.setZero();
  f_s_e.setZero();

  for (auto i = 0; i < GP.getTotalGP(); ++i) {
    prec jacobian = elem->getBeamJacobian(pointers, GP.getXi(i));
    elem->getH1BeamShapes(pointers, this->intOrderDisp, jacobian, shape,
                          shapeDerivative, GP.getXi(i));

    auto &vert1 = elem->getVertex(pointers, 0);
    auto &vert2 = elem->getVertex(pointers, 3);
    Types::Matrix3X<prec> coors;
    coors.resize(3, 2);
    coors.block(0, 0, 3, 1) = vert1.getCoordinates();
    coors.block(0, 1, 3, 1) = vert2.getCoordinates();

    BMatrix = this->getNonBMatrix(shapeDerivative, shape, R01, R02, coors,
                                  solution, numDofs, numberBeamNodes);
    strains_g = this->getstrain_Vector(shapeDerivative, shape, R01, R02, coors,
                                       solution, numberBeamNodes);
    strains = epsilon_hat;
    sigma = sigma_hat;

    K_sigmaMatrix =
        this->getK_sigmaMatrix(shapeDerivative, shape, R01, R02, coors,
                               solution, sigma, numDofs, numberBeamNodes);

    Kg += K_sigmaMatrix * jacobian * GP.getWeight(i);
    G += Amat.transpose() * BMatrix * jacobian * GP.getWeight(i);
    // G += BMatrix.transpose() * Amat * jacobian * GP.getWeight(i);
    F -= Amat.transpose() * Amat * jacobian * GP.getWeight(i);
    H += Amat.transpose() * dmat * Amat * jacobian * GP.getWeight(i);
    f_e_e += Amat.transpose() * dmat * strains * jacobian * GP.getWeight(i);
    f_s_e += Amat.transpose() * strains_g * jacobian * GP.getWeight(i);
  }

  f_i = G.transpose() * sigma_hat;
  f_e = f_e_e + F * sigma_hat;
  f_s = f_s_e + F.transpose() * epsilon_hat;

  // stiffness = Kg +  G * (H.transpose()).inverse() * L * H.inverse()*
  // G.transpose();
  stiffness =
      Kg + G.transpose() * F.inverse() * H * (F.inverse()).transpose() * G;

  residual = G.transpose() *
             (sigma_hat + F.inverse() * H * (F.inverse()).transpose() * f_s -
              F.inverse() * f_e);

  G_previous = G;
  F_previous = F;
  H_previous = H;
  f_s_previous = f_s;
  f_e_previous = f_e;
};

void EL104_TimoshenkoPrism::toParaviewAdaper(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem,
    vtkPlotInterface &paraviewAdapter, ParaviewSwitch control) {
  switch (elem->getType()) {
  case FiniteElement::Elementtypes::LinearPrism: {
    FiniteElement::LinearPrism *passElem;
    passElem = static_cast<FiniteElement::LinearPrism *>(elem);
    this->toParaviewPrism(pointers, passElem, paraviewAdapter, control);
  } break;
  default:
    throw std::runtime_error("Paraview not implemented in this "
                             "elementformulation for the given element type.");
  }
}

auto EL104_TimoshenkoPrism::getHistoryDataStructure()
    -> const HistoryDataStructure & {
  return m_historyDataStructure;
}

auto EL104_TimoshenkoPrism::getNumberOfIntergrationPoints(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> indexType {

  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(this->intOrderDisp);
  return GP.getTotalGP();
}

void EL104_TimoshenkoPrism::toParaviewPrism(PointerCollection &pointers,
                                            FiniteElement::LinearPrism *elem,
                                            vtkPlotInterface &paraviewAdapter,
                                            ParaviewSwitch control) {
  switch (control) {
  case ParaviewSwitch::Mesh: {
    indexType matNum = elem->getMaterial()->getNumber();
    Types::Vector3<prec> coors(3);
    std::vector<indexType> cell(6);
    for (indexType i = 0; i < 6; ++i) {
      auto &Vert = elem->getVertex(pointers, i);
      coors = Vert.getCoordinates();
      paraviewAdapter.addPoint(0, matNum, Vert.getId(), coors(0), coors(1),
                               coors(2));
      cell[i] = Vert.getId();
    }
    int celltype = VTK_WEDGE;
    indexType dummy = 1;
    indexType nn = 6;
    paraviewAdapter.addCell(0, matNum, elem->getId(), dummy, cell, nn,
                            celltype);

  } break;
  case ParaviewSwitch::Solution: {
    indexType matNum = elem->getMaterial()->getNumber();
    std::vector<Types::Matrix33<prec>> bases;
    bases.push_back(elem->getStartTriad(pointers));
    bases.push_back(elem->getEndTriad(pointers));

    Geometry::Vertex *V1, *V2;
    V1 = &elem->getVertex(pointers, 0);
    V2 = &elem->getVertex(pointers, 3);

    std::vector<Types::Vector3<prec>> positions;
    positions.push_back(V1->getCoordinates());
    positions.push_back(V2->getCoordinates());

    std::vector<DegreeOfFreedom *> dispDofs, rotDofs;

    elem->getH1BeamDofs(pointers, dispDofs, this->meshIdDisp,
                        this->intOrderDisp);
    elem->getH1BeamDofs(pointers, rotDofs, this->meshIdRot, this->intOrderDisp);

    Types::VectorX<prec> rotSolution, dispSolution;
    elem->getSolution(pointers, dispDofs, dispSolution);
    elem->getSolution(pointers, rotDofs, rotSolution);

    Types::Matrix33<prec> skewMat;

    indexType vertNumber = 0;
    indexType localVertNumber = 0;
    std::vector<prec> vsol, vertRot;
    Types::Vector3<prec> vertSol;

    for (auto i = 0; i < 2; ++i) {
      Types::Vector3<prec> omega;
      omega = rotSolution.block(3 * localVertNumber, 0, 3, 1);
      vertRot.clear();
      vertRot.push_back(omega(0));
      vertRot.push_back(omega(1));
      vertRot.push_back(omega(2));

      Types::Matrix33<prec> triads;
      triads = bases[i];

      // Vertex 1
      vsol.clear();
      vsol.push_back(dispSolution(0 + 3 * localVertNumber));
      vsol.push_back(dispSolution(1 + 3 * localVertNumber));
      vsol.push_back(dispSolution(2 + 3 * localVertNumber));

      Geometry::Vertex *Vert;
      Vert = &elem->getVertex(pointers, vertNumber);

      paraviewAdapter.setPointData(0, matNum, Vert->getId(), vsol, 3,
                                   paraviewNames::DisplacementName());
      paraviewAdapter.setPointData(0, matNum, Vert->getId(), vertRot, 3,
                                   paraviewNames::RotationName());

      // Vert 2
      skewMat = Math::Geometry::skewMatrix(omega);
      vertSol = dispSolution.block(3 * localVertNumber, 0, 3, 1);

      vertSol += skewMat * bases[localVertNumber].block(0, 1, 3, 1);

      Geometry::Vertex *Vert2;
      Vert2 = &elem->getVertex(pointers, vertNumber + 1);

      vsol.clear();
      vsol.push_back(vertSol(0));
      vsol.push_back(vertSol(1));
      vsol.push_back(vertSol(2));

      Types::Vector3<prec> axis2;
      axis2 = triads.block(0, 1, 3, 1);
      Types::Vector3<prec> icoor, ccoor;
      icoor = Vert->getCoordinates() + axis2;

      auto R = Math::Geometry::getRotationMatrix(omega);
      axis2 = R * axis2;
      ccoor =
          Vert->getCoordinates() + dispSolution.block(3 * i, 0, 3, 1) + axis2;

      Types::Vector3<prec> diffdisp = ccoor - icoor;
      std::vector<prec> tempDisp;
      tempDisp.clear();
      for (auto kk = 0; kk < 3; ++kk) {
        tempDisp.push_back(diffdisp(kk));
      }

      paraviewAdapter.setPointData(0, matNum, Vert2->getId(), tempDisp, 3,
                                   paraviewNames::DisplacementName());
      paraviewAdapter.setPointData(0, matNum, Vert2->getId(), vertRot, 3,
                                   paraviewNames::RotationName());
      // Vert 3

      vertSol = dispSolution.block(3 * localVertNumber, 0, 3, 1);

      vertSol += skewMat * bases[localVertNumber].block(0, 2, 3, 1);
      Vert2 = &elem->getVertex(pointers, vertNumber + 2);

      vsol.clear();
      vsol.push_back(vertSol(0));
      vsol.push_back(vertSol(1));
      vsol.push_back(vertSol(2));

      axis2 = triads.block(0, 2, 3, 1);
      icoor = Vert->getCoordinates() + axis2;
      axis2 = R * axis2;
      ccoor =
          Vert->getCoordinates() + dispSolution.block(3 * i, 0, 3, 1) + axis2;

      diffdisp = ccoor - icoor;
      tempDisp.clear();
      for (auto kk = 0; kk < 3; ++kk) {
        tempDisp.push_back(diffdisp(kk));
      }

      paraviewAdapter.setPointData(0, matNum, Vert2->getId(), tempDisp, 3,
                                   paraviewNames::DisplacementName());
      paraviewAdapter.setPointData(0, matNum, Vert2->getId(), vertRot, 3,
                                   paraviewNames::RotationName());

      vertNumber += 3;
      localVertNumber++;
    }

  } break;
  case HierAMuS::ParaviewSwitch::ProjectedValues: {
    indexType matNum = elem->getMaterial()->getNumber();

    std::vector<DegreeOfFreedom *> Dofs;
    Types::MatrixXX<prec> stiffness;
    Types::VectorX<prec> residual;

    this->setTangentResidual(pointers, elem, stiffness, residual, Dofs);

    indexType vertNumber = 0;
    for (auto i = 0; i < 2; ++i) {

      auto &Vert = elem->getVertex(pointers, vertNumber);
      vertNumber += 3;

      std::vector<prec> transReact, rotReact;
      for (auto j = 0; j < 3; ++j) {
        transReact.push_back(residual(j + 3 * i));
        rotReact.push_back(residual(j + 3 * i + 6));
      }
      paraviewAdapter.SumPointData(0, matNum, transReact, Vert.getId(), 3,
                                   paraviewNames::ReactionsTransName());
      paraviewAdapter.SumPointData(0, matNum, rotReact, Vert.getId(), 3,
                                   paraviewNames::ReactionsRotName());
    }
  } break;
  default:
    break;
  }
}

Types::Matrix6X<prec>
EL104_TimoshenkoPrism::getLocalStrainStressInterpolation() {

  Types::Matrix6X<prec> Amat;
  Amat.resize(6, 6);
  Amat.setZero();

  for (auto i = 0; i < 6; ++i) {
    Amat(i, i) = prec(1);
  }

  return Amat;
}

const HistoryDataStructure
EL104_TimoshenkoPrism::m_historyDataStructure({{6,12},{6,6},{6,6},{6,1},{6,1},{6,1},{6,1}},{});

} // namespace HierAMuS::Elementformulations
