// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include <elementFormulations/EL307_VolumeConstraint.h>
#include <pointercollection/pointercollection.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include "MatrixTypes.h"
#include "control/ParameterList.h"

#include "finiteElements/GenericFiniteElement.h"
#include "finiteElements/VolumeConstraint.h"
#include "materials/GenericMaterialFormulation.h"

#include "math/MatrixOperations.h"



namespace HierAMuS::Elementformulations {


EL307_VolumeConstraint::EL307_VolumeConstraint(PointerCollection *ptrCol)
    : GenericElementFormulation(ptrCol) {}

EL307_VolumeConstraint::~EL307_VolumeConstraint() {}


void EL307_VolumeConstraint::readData(PointerCollection &pointers,
                                    ParameterList &list) {

  m_mode = list.getIndexVal("mode");
  m_meshIdDisp = list.getIndexVal("meshiddisp");
  m_meshIdLam = list.getIndexVal("meshIdLam");
  m_intOrderDisp = list.getIndexVal("disporder");
  m_center = list.getIndexVal("center");

  m_K = list.getPrecVal("stiffness");


  auto Log = pointers.getSPDLogger();

  Log.info("\n{:-<100}\n"
                "*   Element 307, specified Options\n"
                "    Mesh id for displacement nodes:         {:>12}\n"
                "    Shape function order for displacements: {:>12}\n"
                "{:-<100}\n",
                "",
                this->m_meshIdDisp,
                this->m_intOrderDisp,
                "");


  this->messageUnprocessed(pointers, list, "EL307_VolumeConstraint");
}

void EL307_VolumeConstraint::setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
{
  FiniteElement::VolumeConstraint *volelem =
      dynamic_cast<FiniteElement::VolumeConstraint *>(elem);

  volelem->setH1Shapes(pointers, m_meshIdDisp, m_intOrderDisp);
  volelem->setVertexNodes(pointers, m_meshIdLam);


  Types::Vector3<prec> coor = pointers.getGeometryData()->getxMax();
  m_xmax = coor(0);
}

void EL307_VolumeConstraint::AdditionalOperations(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) {

  FiniteElement::VolumeConstraint *volelem =
      dynamic_cast<FiniteElement::VolumeConstraint *>(elem);

  if (m_mode == 4) { // Thermal boundary conditions
    auto Nodes = volelem->getVertexNodes(pointers, m_meshIdLam);
    Nodes[0]->getDegreeOfFreedom(2).setStatus(dofStatus::inactive);
  }
}

auto EL307_VolumeConstraint::getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
-> std::vector<DegreeOfFreedom *>  {

  FiniteElement::VolumeConstraint *volelem =
      dynamic_cast<FiniteElement::VolumeConstraint *>(elem);

  std::vector<DegreeOfFreedom *> Dofs;
  volelem->getH1Dofs(pointers, Dofs, m_meshIdDisp, m_intOrderDisp);
  volelem->getVertexDofs(pointers, Dofs, m_meshIdLam);
  return Dofs;
  ;
}

void EL307_VolumeConstraint::setTangentResidual(PointerCollection& pointers,
                                              FiniteElement::GenericFiniteElement *elem,
                                              Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
                                              Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs)
{
  FiniteElement::VolumeConstraint *volElement =
      dynamic_cast<FiniteElement::VolumeConstraint *>(elem);

  switch(this->m_mode){
    case 1:
    this->setTangentResidual_XDir(pointers, volElement, stiffness,
                                  residual, Dofs);
    break;
  case 4:
    this->setTangentResidual_TempGradient(pointers, volElement, stiffness,
                                          residual, Dofs);
    break;
  }


}

auto EL307_VolumeConstraint::getNumberOfIntergrationPoints(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> indexType {
  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(this->m_intOrderDisp);
  return GP.getTotalGP();
}

void EL307_VolumeConstraint::toParaviewAdaper(PointerCollection &pointers,
                                FiniteElement::GenericFiniteElement *elem,
                                vtkPlotInterface &paraviewAdapter,
                                ParaviewSwitch control)
{

}

void EL307_VolumeConstraint::setTangentResidual_XDir(
  PointerCollection &pointers,
  FiniteElement::VolumeConstraint* volElement,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs)
{
  Dofs.clear();
  volElement->getH1Dofs(pointers, Dofs, m_meshIdDisp, m_intOrderDisp);
  volElement->getVertexDofs(pointers, Dofs, m_meshIdLam);

  indexType numDofs = Dofs.size();
  stiffness.resize(numDofs, numDofs);
  residual.resize(numDofs);
  stiffness.setZero();
  residual.setZero();

  Types::VectorX<prec> solution = volElement->getSolution(pointers, Dofs);

  auto GP = volElement->getIntegrationPoints(pointers);
  GP.setOrder(m_intOrderDisp * m_intOrderDisp);
  Types::VectorX<prec> Neps;
  Neps.resize(numDofs-3);
  Neps.setZero();

  Types::Vector3<prec> p;
  indexType numDispNodes = numDofs / 3 - 1;
  for (auto gp:GP)
  {
    Types::MatrixXX<prec> jaco = volElement->getJacobian(pointers, gp);
    auto shapes = volElement->getH1Shapes(pointers, m_intOrderDisp, jaco, gp);
    Types::Vector3<prec> coor = volElement->getVolumeCoordinates(pointers, gp);
    prec dV = jaco.determinant() * gp.weight;
    prec NLam;
    if (m_center <= 0.1)
    {
      NLam = coor(0);
    } else
    {
      coor(0) > 0 ? NLam = coor(0) - m_xmax : NLam = coor(0) + m_xmax;
    }
    for (indexType i=0;i<numDispNodes;++i)
    {
      Neps(3 * i) = shapes.shapeDeriv(0,i);
    }
    //std::cout << "Nlam: " << NLam << "  coor: " << coor.transpose()
    //          << std::endl;
    p(0) = NLam;
    p(1) = NLam * coor(1);
    p(2) = NLam * coor(2);
    stiffness.block(0, numDofs - 3, numDofs - 3, 3) +=
        Neps * p.transpose() * dV;
    //std::cout << "Stiffness part constraint:\n" << Neps * p.transpose() * dV << std::endl;
    //stiffness.block( numDofs - 3, 0, 3, numDofs - 3) += p * Neps.transpose() * dV;

  }
  stiffness += stiffness.transpose().eval();
  stiffness *= m_K;
  residual += stiffness * solution;
}

void EL307_VolumeConstraint::setTangentResidual_TempGradient(
    PointerCollection &pointers, FiniteElement::VolumeConstraint *volElement,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs)
{

  Dofs.clear();
  volElement->getH1Dofs(pointers, Dofs, m_meshIdDisp, m_intOrderDisp);
  volElement->getVertexDofs(pointers, Dofs, m_meshIdLam);

  
  indexType numDofs = Dofs.size();
  stiffness.resize(numDofs, numDofs);
  residual.resize(numDofs);
  stiffness.setZero();
  residual.setZero();

  Types::VectorX<prec> solution = volElement->getSolution(pointers, Dofs);

  auto GP = volElement->getIntegrationPoints(pointers);
  GP.setOrder(m_intOrderDisp * m_intOrderDisp);
  Types::Matrix3X<prec> Neps;
  Neps.resize(3,numDofs - 3);
  Neps.setZero();
  
  indexType numDispNodes = numDofs / 3 - 1;
  for (auto gp : GP) {
    Types::MatrixXX<prec> jaco = volElement->getJacobian(pointers, gp);
    auto shapes = volElement->getH1Shapes(pointers, m_intOrderDisp, jaco, gp);
    Types::Vector3<prec> coor = volElement->getVolumeCoordinates(pointers, gp);
    prec dV = jaco.determinant() * gp.weight;
    prec NLam;
    if (m_center <= 0.1) {
      NLam = coor(0);
    } else {
      coor(0) > 0 ? NLam = coor(0) - m_xmax : NLam = coor(0) + m_xmax;
    }
    for (indexType i = 0; i < numDispNodes; ++i) {
      Neps(0,3 * i) = shapes.shapeDeriv(1, i);
      Neps(1,3 * i) = shapes.shapeDeriv(2, i);
    }
    // std::cout << "Nlam: " << NLam << "  coor: " << coor.transpose()
    //           << std::endl;
    dV *= NLam;
    stiffness.block(0, numDofs - 3, numDofs - 3, 3) +=
        Neps.transpose() * dV;
    //std::cout << "Stiffness part constraint:\n" << Neps * p.transpose() * dV
    // << std::endl; stiffness.block( numDofs - 3, 0, 3, numDofs - 3) += p *
    // Neps.transpose() * dV;
  }
  stiffness += stiffness.transpose().eval();
  stiffness *= m_K;
  residual += stiffness * solution;
  //std::cout << stiffness << std::endl;

}

auto EL307_VolumeConstraint::getBMatrixLinear(PointerCollection& pointers, Geometry::H1Shapes &shapes, IntegrationPoint &intPoint, FiniteElement::beamInterfaceElement3D &elem) -> Types::Matrix6X<prec>
{
  indexType numShapes = shapes.shapes.rows();

  Types::Vector3<prec> localCoor = elem.getLocalCoordinate(pointers, intPoint);

  Types::Matrix6X<prec> BMat;
  BMat.resize(6, numShapes*3);
  BMat.setZero();

  Types::Matrix33<prec> R0T = elem.getRotationR0().transpose();
  Types::Vector3T<prec> A1 = R0T.block(0, 0, 1, 3);
  Types::Vector3T<prec> A2 = R0T.block(1, 0, 1, 3);
  Types::Vector3T<prec> A3 = R0T.block(2, 0, 1, 3);

  for (auto i=0;i<numShapes-1;++i){
    BMat.block(0,i*3,1,3) = shapes.shapeDeriv(0,i)*A1;
    BMat.block(1,i*3,1,3) = shapes.shapeDeriv(1,i)*A2;
    BMat.block(2,i*3,1,3) = shapes.shapeDeriv(2,i)*A3;

    BMat.block(3,i*3,1,3) = shapes.shapeDeriv(1,i)*A1;
    BMat.block(3,i*3,1,3) += shapes.shapeDeriv(0,i)*A2;

    BMat.block(4,i*3,1,3) = shapes.shapeDeriv(2,i)*A1;
    BMat.block(4,i*3,1,3) += shapes.shapeDeriv(0,i)*A3;

    BMat.block(5,i*3,1,3) = shapes.shapeDeriv(2,i)*A2;
    BMat.block(5,i*3,1,3) += shapes.shapeDeriv(1,i)*A3;

  }

  indexType pos = numShapes-1;
  pos *= 3;

  BMat.block(0, pos, 1, 3) = shapes.shapeDeriv(0,numShapes-1)*(localCoor(2)*A2-localCoor(1)*A3);

  BMat.block(3, pos, 1, 3) = -shapes.shapeDeriv(0,numShapes-1)*localCoor(2)*A1;
  BMat.block(3, pos, 1, 3) += -shapes.shapes(numShapes-1)*A3;

  BMat.block(4, pos, 1, 3) = shapes.shapeDeriv(0,numShapes-1)*localCoor(1)*A1;
  BMat.block(4, pos, 1, 3) += shapes.shapes(numShapes-1)*A2;

  //std::cout << "R0\n" << R0T << std::endl;

  return BMat;
}



} // namespace HierAMuS::Elementformulations
