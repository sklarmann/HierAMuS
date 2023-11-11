// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include <elementFormulations/EL302_BeamCoupling3D.h>
#include <pointercollection/pointercollection.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include "MatrixTypes.h"
#include "control/ParameterList.h"

#include "datatypes.h"
#include "finiteElements/GenericFiniteElement.h"
#include "finiteElements/beamInterfaceElement3D.h"
#include "materials/GenericMaterialFormulation.h"
#include "materials/Material.h"

#include "math/MatrixOperations.h"

#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/bundled/ranges.h"

namespace HierAMuS::Elementformulations {


EL302_BeamCoupling3D::EL302_BeamCoupling3D(PointerCollection *ptrCol)
    : GenericElementFormulationInterface(ptrCol) {}

EL302_BeamCoupling3D::~EL302_BeamCoupling3D() {}


void EL302_BeamCoupling3D::readData(PointerCollection &pointers,
                                    ParameterList &list) {
  m_meshIdDisp = list.getIndexVal("meshiddisp");
  m_meshIdRot = list.getIndexVal("meshidrot");
  m_intOrderDisp = list.getIndexVal("disporder");
  m_mode = list.getIndexVal("mode");
  m_warpType = list.getIndexVal("warptype");

  if (list.hasParameter("warpBounNodes"))
  {
    Types::MatrixXX<indexType> verts = list.getIndexMatrix("warpBounNodes");
    for (indexType i=0;i<3;++i)
    {
      m_warpBounVerts.push_back(verts(i));
    }
  }

  
  auto &Log = pointers.getSPDLogger();

  Log.info("\n{:-<100}\n"
                "*   Element 302, specified Options\n"
                "    Mesh id for displacement nodes:         {:>12}\n"
                "    Shape function order for displacements: {:>12}\n"
                "{:-<100}\n",
                "",
                this->m_meshIdDisp,
                this->m_intOrderDisp,
                "");

  if (m_warpBounVerts.size()!=0)
  {
    Log.info("    Constraining local verts:         {}\n", fmt::join(m_warpBounVerts,","));

    
  }

  this->messageUnprocessed(pointers, list, "EL302_BeamCoupling3D");
}

void EL302_BeamCoupling3D::setDegreesOfFreedom(
    PointerCollection &pointers, FiniteElement::beamInterfaceElement3D &elem) {
  FiniteElement::beamInterfaceElement3D *interfaceElement =
      dynamic_cast<FiniteElement::beamInterfaceElement3D *>(&elem);

  interfaceElement->setH1Shapes(pointers, m_intOrderDisp, m_meshIdDisp,
                                m_meshIdRot);
  interfaceElement->setWarpingType(m_warpType);

}

void EL302_BeamCoupling3D::AdditionalOperations(PointerCollection& pointers, FiniteElement::beamInterfaceElement3D &elem)
{


  if (m_warpBounVerts.size()!=0)
    elem.setWarpDofConstraintVertices(m_warpBounVerts);
  elem.computeWarpingShapes(pointers);

  elem.print(pointers);
  

}

auto EL302_BeamCoupling3D::getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
-> std::vector<DegreeOfFreedom *>  {
  FiniteElement::beamInterfaceElement3D *interfaceElement =
      dynamic_cast<FiniteElement::beamInterfaceElement3D *>(elem);
  return interfaceElement->getH1Dofs(pointers, m_meshIdDisp, m_meshIdRot,
                                     m_intOrderDisp);
  ;
}

void EL302_BeamCoupling3D::setTangentResidual(PointerCollection& pointers, FiniteElement::beamInterfaceElement3D &elem,
                                              Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
                                              Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs)
{
  FiniteElement::beamInterfaceElement3D *interfaceElement = &elem;

  switch(this->m_mode){
    case 1:
    this->setTangentResidual_PureDisp(pointers, interfaceElement, stiffness,
                                      residual, Dofs);
    break;
    case 2:
      this->setTangentResidual_HuWashizu(pointers, interfaceElement, stiffness,
                                         residual, Dofs);
    break;
    case 3:
      this->setTangentResidual_HuWashizu_Warping(pointers, interfaceElement,
                                                 stiffness, residual, Dofs);
    break;
    case 4:
      this->setTangentResidual_Disp_Warping_V1(pointers, interfaceElement,
                                                 stiffness, residual, Dofs);
    break;
  }
  


}

auto EL302_BeamCoupling3D::getNumberOfIntergrationPoints(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> indexType {
  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(this->m_intOrderDisp);
  return GP.getTotalGP();
}

void EL302_BeamCoupling3D::toParaviewAdaper(PointerCollection &pointers,
                                FiniteElement::beamInterfaceElement3D &elem,
                                vtkPlotInterface &paraviewAdapter,
                                ParaviewSwitch control)
{

  int matNum = static_cast<int>(elem.getMaterial()->getNumber());

  switch(control) {
  case ParaviewSwitch::Mesh: {
    elem.geometryToParaview(pointers, paraviewAdapter, 0, matNum);
  } break;
  case ParaviewSwitch::Solution: {
    elem.H1SolutionToParaview(pointers,paraviewAdapter,0,matNum,m_meshIdDisp,m_intOrderDisp,paraviewNames::DisplacementName());
  } break;
  default:
    break;
  }
}

void EL302_BeamCoupling3D::setTangentResidual_PureDisp(PointerCollection& pointers,
                                                       FiniteElement::beamInterfaceElement3D *interfaceElement,
                                                       Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
                                                       Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs)
{
  Dofs = interfaceElement->getH1Dofs(pointers, m_meshIdDisp, m_meshIdRot,
                                     m_intOrderDisp);
  stiffness.resize(Dofs.size(), Dofs.size());
  residual.resize(Dofs.size());
  stiffness.setZero();
  residual.setZero();


  Types::VectorX<prec> solution = interfaceElement->getSolution(pointers, Dofs);

  Materials::MaterialTransferData materialTransfer;
  materialTransfer.strains.resize(6);

  auto GP = interfaceElement->getIntegrationPoints(pointers);
  GP.setOrder(m_intOrderDisp*2);

  Types::VectorX<prec> templsg (Dofs.size());
  //templsg.setZero();
  //templsg(Dofs.size()-6) = prec(1);
  //templsg(Dofs.size()-5) = prec(1);
  //templsg(Dofs.size()-4) = prec(1);

  for(auto i:GP){
    auto jacobi = interfaceElement->getJacobian(pointers, i);
    auto dispShapes = interfaceElement->getH1Shapes(pointers, m_intOrderDisp,
                                                    m_meshIdDisp, jacobi, i);

    prec dV = jacobi.determinant()*i.weight;

    Types::Matrix6X<prec> BMat =
        this->getBMatrixLinear(pointers, dispShapes, i, *interfaceElement);
    //auto testeps = BMat*templsg;
    //std::cout << "   Test strains: " << testeps.transpose() << std::endl;

    materialTransfer.strains = BMat*solution;
    interfaceElement->getMaterialFormulation(pointers, i)->getMaterialData(pointers, materialTransfer,i);

    stiffness += BMat.transpose()*materialTransfer.materialTangent*BMat*dV;
    residual += BMat.transpose()*materialTransfer.stresses*dV;

  }
  //#pragma omp critical
  //{
  //  Types::VectorX<prec> evals = stiffness.eigenvalues().real();
  //  std::sort(evals.data(), evals.data()+evals.size()); 
  //  std::cout << "  Eigenvalues of stiffness matrix: " << evals << std::endl;
  //}
}

void EL302_BeamCoupling3D::setTangentResidual_HuWashizu(PointerCollection& pointers,
                                                        FiniteElement::beamInterfaceElement3D *interfaceElement,
                                                        Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
                                                        Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs)
{
  Dofs = interfaceElement->getH1Dofs(pointers, m_meshIdDisp, m_meshIdRot,
                                     m_intOrderDisp);
  stiffness.resize(Dofs.size(), Dofs.size());
  residual.resize(Dofs.size());
  stiffness.setZero();
  residual.setZero();

  indexType nEps = interfaceElement->getNumberOfInnerDofs(this->m_intOrderDisp);

  Types::VectorX<prec> solution = interfaceElement->getSolution(pointers, Dofs);

  Materials::MaterialTransferData materialTransfer;
  materialTransfer.strains.resize(6);

  Types::MatrixXX<prec> H, G, L;
  H.resize(nEps, nEps);
  G.resize(nEps, nEps);
  L.resize(Dofs.size(), nEps);
  H.setZero();
  G.setZero();
  L.setZero();

  auto GP = interfaceElement->getIntegrationPoints(pointers);
  GP.setOrder(m_intOrderDisp*2);

  for(auto i:GP){
    auto jacobi = interfaceElement->getJacobian(pointers, i);
    auto dispShapes = interfaceElement->getH1Shapes(pointers, m_intOrderDisp,
                                                    m_meshIdDisp, jacobi, i);

    prec dV = jacobi.determinant()*i.weight;

    auto A = interfaceElement->getL2ShapesStressStrain(
        pointers, i, this->m_intOrderDisp, jacobi);

    Types::Matrix6X<prec> BMat =
        this->getBMatrixLinear(pointers, dispShapes, i, *interfaceElement);

    //std::cout << A << "\n\n" << std::endl;

    materialTransfer.strains = BMat*solution;
    interfaceElement->getMaterialFormulation(pointers, i)->getMaterialData(pointers, materialTransfer,i);

    // Types::MatrixXX<prec> Aq = A.shapes.transpose()*A.shapes*dV;
    // std::cout << i.sectionNumber << std::endl;
    // std::cout << Aq << std::endl;

    G.block(A.pos, A.pos, A.numShapes, A.numShapes)+= A.shapes.transpose()*A.shapes*dV;
    L.block(0,A.pos,L.rows(), A.numShapes) += BMat.transpose()*A.shapes*dV;
    H.block(A.pos, A.pos, A.numShapes, A.numShapes)+= A.shapes.transpose()*materialTransfer.materialTangent*A.shapes*dV;


  }


  G = G.inverse();
  Types::MatrixXX<prec> TempMat = G*H*G.transpose();

  stiffness = L*TempMat*L.transpose();
  residual = stiffness*solution;

    //stiffness += BMat.transpose()*materialTransfer.materialTangent*BMat*dV;
    //residual += BMat.transpose()*materialTransfer.stresses*dV;
}

void EL302_BeamCoupling3D::setTangentResidual_HuWashizu_Warping(PointerCollection& pointers,
                                                                FiniteElement::beamInterfaceElement3D *interfaceElement,
                                                                Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
                                                                Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs)
{
  Dofs = interfaceElement->getH1Dofs(pointers, m_meshIdDisp, m_meshIdRot,
                                     m_intOrderDisp);
  stiffness.resize(Dofs.size(), Dofs.size());
  residual.resize(Dofs.size());
  stiffness.setZero();
  residual.setZero();

  indexType nEps = interfaceElement->getNumberOfInnerDofs(this->m_intOrderDisp);
  indexType nwarp = interfaceElement->getNumberOfInnerWarpingDofs();

  Types::VectorX<prec> solution = interfaceElement->getSolution(pointers, Dofs);

  Materials::MaterialTransferData materialTransfer;
  materialTransfer.strains.resize(6);

  Types::MatrixXX<prec> H, G, L, M;
  H.resize(nEps, nEps);
  M.resize(nEps, nEps);
  G.resize(Dofs.size(), nEps);
  L.resize(nEps, nwarp);
  H.setZero();
  G.setZero();
  L.setZero();
  M.setZero();

  auto GP = interfaceElement->getIntegrationPoints(pointers);
  GP.setOrder(m_intOrderDisp*2);

  for(auto i:GP){
    auto jacobi = interfaceElement->getJacobian(pointers, i);
    auto dispShapes = interfaceElement->getH1Shapes(pointers, m_intOrderDisp,
                                                    m_meshIdDisp, jacobi, i);

    prec dV = jacobi.determinant()*i.weight;

    auto A = interfaceElement->getL2ShapesStressStrain(
        pointers, i, this->m_intOrderDisp, jacobi);

    Types::Matrix6X<prec> BMat =
        this->getBMatrixLinear(pointers, dispShapes, i, *interfaceElement);

    auto warpMat = interfaceElement->getWarpingMatrix(pointers, i, m_intOrderDisp, jacobi);


  //std::cout << warpMat << "\n" << std::endl;


    //std::cout << A << "\n\n" << std::endl;

    materialTransfer.strains = BMat*solution;
    interfaceElement->getMaterialFormulation(pointers, i)->getMaterialData(pointers, materialTransfer,i);

    H.block(A.pos, A.pos, A.numShapes, A.numShapes)+= A.shapes.transpose()*A.shapes*dV;
    G.block(0,A.pos,BMat.cols(), A.numShapes) += BMat.transpose()*A.shapes*dV;
    M.block(A.pos, A.pos, A.numShapes, A.numShapes)+= A.shapes.transpose()*materialTransfer.materialTangent*materialTransfer.materialTangent*A.shapes*dV;
    L.block(A.pos, 0, A.numShapes, nwarp) +=
        A.shapes.transpose() * warpMat * dV;

    //std::cout << "A:\n" << A.shapes << "\n\nWarp: \n" << warpMat << "\n\n" << std::endl;


  }

  //std::cout << H.eigenvalues() << std::endl;
  //Types::MatrixXX<prec> MinvH=Math::AinvTimesB(M, H);
  //std::cout << H << std::endl;
  {
    Eigen::BDCSVD<Types::MatrixXX<prec>> lu_decomp(H);
    auto rank = lu_decomp.rank();
    std::cout << "Rank of H: " << rank << std::endl;
    std::cout << "Cols of H: " << H.cols() << std::endl;
  }

  H = H.inverse();
  for (auto i = 0; i < H.rows(); ++i) {
    for (auto j = i; j < H.cols(); ++j) {
      prec tt = (H(i, j) + H(j, i)) * prec(0.5);
      H(i, j) = tt;
      H(j, i) = tt;
    }
  }
  //H = H.transpose() * M * H;

  Types::MatrixXX<prec> Q = H.transpose() * M * H;
  //Eigen::LDLT<Types::MatrixXX<prec>> HLDL;
  //HLDL.compute(H);
  //= H.selfadjointView<Eigen::Upper>().ldlt();

  //Types::MatrixXX<prec> Q = HLDL.solve(HLDL.solve(M).transpose());

  for (auto i = 0; i < Q.rows(); ++i) {
    for (auto j = i; j < Q.cols(); ++j) {
      prec tt = (Q(i, j) + Q(j, i)) * prec(0.5);
      Q(i, j) = tt;
      Q(j, i) = tt;
    }
  }
  {
    Eigen::BDCSVD<Types::MatrixXX<prec>> lu_decomp(Q);
    auto rank = lu_decomp.rank();
    std::cout << "Rank of Q: " << rank << std::endl;
    std::cout << "Cols of Q: " << Q.cols() << std::endl;
  }


  //Types::MatrixXX<prec> QinvL = Math::AinvTimesB(H, L);

  Types::MatrixXX<prec> tempMat = G * Q* L;
  Types::MatrixXX<prec> tempMat2 =
      L.transpose() * Q * L;
  {
    Eigen::BDCSVD<Types::MatrixXX<prec>> lu_decomp(tempMat2);
    auto rank = lu_decomp.rank();
    std::cout << "Rank of tempMat2: " << rank << std::endl;
    std::cout << "Cols of tempMat2: " << tempMat2.cols() << std::endl;
  }
  //Types::MatrixXX<prec> tempMat = G * QinvL;
  //Types::MatrixXX<prec> tempMat2 = L.transpose() * QinvL;
  //std::cout << tempMat2.eigenvalues() << std::endl;

  //Types::MatrixXX<prec> GT = G.transpose();
  //Types::MatrixXX<prec> QinvGT = Math::AinvTimesB(H, GT);

  stiffness = G * Q * G.transpose() -
              tempMat * tempMat2.inverse() * tempMat.transpose();
  //(tempMat *tempMat2.selfadjointView<Eigen::Upper>().ldlt().solve(
  //                           tempMat.transpose()));
  for (auto i = 0; i < stiffness.rows(); ++i) {
    for (auto j = i; j < stiffness.cols(); ++j) {
      prec tt = (stiffness(i, j) + stiffness(j, i)) * prec(0.5);
      stiffness(i, j) = tt;
      stiffness(j, i) = tt;
    }
  }
  {
    Eigen::BDCSVD<Types::MatrixXX<prec>> lu_decomp(stiffness);
    auto rank = lu_decomp.rank();
    std::cout << "Rank of stiffness: " << rank << std::endl;
    std::cout << "Cols of stiffness: " << stiffness.cols() << std::endl;
  }

  residual = stiffness.selfadjointView<Eigen::Upper>() * solution;


#pragma omp critical
  {
  //std::cout << stiffness.eigenvalues() << std::endl;
}
    //stiffness += BMat.transpose()*materialTransfer.materialTangent*BMat*dV;
    //residual += BMat.transpose()*materialTransfer.stresses*dV;
}

void EL302_BeamCoupling3D::setTangentResidual_Disp_Warping_V1(PointerCollection& pointers,
    FiniteElement::beamInterfaceElement3D *interfaceElement,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs)
{
  Dofs = interfaceElement->getH1Dofs(pointers, m_meshIdDisp, m_meshIdRot,
                                     m_intOrderDisp);
  stiffness.resize(Dofs.size(), Dofs.size());
  residual.resize(Dofs.size());
  stiffness.setZero();
  residual.setZero();

  indexType sizeA = Dofs.size();
  indexType sizeB = sizeA - 6;
  Types::MatrixXX<prec> L(sizeB, sizeA);
  Types::MatrixXX<prec> M(sizeB, sizeB);
  L.setZero();
  M.setZero();

  Types::VectorX<prec> solution = interfaceElement->getSolution(pointers, Dofs);

  Materials::MaterialTransferData materialTransfer;
  materialTransfer.strains.resize(6);

  auto GP = interfaceElement->getIntegrationPoints(pointers);
  GP.setOrder(m_intOrderDisp*2);
  //GP.setNumGpPerDir(1, m_intOrderDisp + 1, m_intOrderDisp + 1);


  for(auto i:GP){
    auto jacobi = interfaceElement->getJacobian(pointers, i);
    auto dispShapes = interfaceElement->getH1Shapes(pointers, m_intOrderDisp,
                                                    m_meshIdDisp, jacobi, i);

    prec dV = jacobi.determinant()*i.weight;

    //Types::Matrix6X<prec> BMat =
    //    this->getBMatrixLinear(pointers, dispShapes, i, *interfaceElement);
    
    Types::Matrix6X<prec> BMat = interfaceElement->getBMatrixWarpingLocal2(
        pointers, i, m_meshIdDisp, m_intOrderDisp, jacobi);
    materialTransfer.strains = BMat*solution;
    interfaceElement->getMaterialFormulation(pointers, i)->getMaterialData(pointers, materialTransfer,i);

    //auto BWarp = interfaceElement->getBMatrixWarpingLocal(pointers, i, this->m_meshIdDisp, this->m_intOrderDisp, jacobi);
    //std::cout << "BMatWarp \n" << BWarp << "\n" << std::endl;

    //L += BWarp.transpose()*materialTransfer.materialTangent*BMat*dV;
    //M += BWarp.transpose()*materialTransfer.materialTangent*BWarp*dV;

    stiffness += BMat.transpose()*materialTransfer.materialTangent*BMat*dV;
    //residual += BMat.transpose()*materialTransfer.stresses*dV;

  }

  

  //std::cout << "M ev:\n" << M.eigenvalues() << std::endl;



  //std::cout << "stiffness before ev:\n" << stiffness.eigenvalues() << std::endl;
  //auto temp = (M.lu().solve(L)).eval();
  //stiffness -= L.transpose()*temp;
  residual = stiffness * solution;

  //std::cout << "stiffness after ev:\n" << stiffness.eigenvalues() << std::endl;
  //auto warpSol = -temp*solution;

  //Eigen::IOFormat CleanFmt(6, 0, "; ", "\n", "", ";");

  //for (auto i:GP){
  //  auto jacobi = interfaceElement->getJacobian(pointers, i);
  //  auto dispShapes = interfaceElement->getH1Shapes(pointers, m_intOrderDisp,
  //                                                  m_meshIdDisp, jacobi, i);

  //  prec dV = jacobi.determinant()*i.weight;

  //  Types::Matrix6X<prec> BMat =
  //      this->getBMatrixLinear(pointers, dispShapes, i, *interfaceElement);
  //  

  //  auto BWarp = interfaceElement->getBMatrixWarpingLocal(pointers, i, this->m_meshIdDisp, this->m_intOrderDisp, jacobi);
  //  Types::Matrix6X<prec> BMatWarp2;
  //  BMatWarp2.resize(6, BWarp.cols()-6);
  //  indexType cc = 0;
  //  for (indexType j=0;j<BWarp.cols();++j)
  //  {
  //    if (removeCols.find(j)==removeCols.end()) 
  //    {
  //      //BMatWarp2.block<6,1>(0,cc) = BWarp.block<6,1>(0,j);
  //      BMatWarp2(Eigen::all,cc) = BWarp(Eigen::all,j);
  //      ++cc;
  //    }
  //  }
  //  
  //  materialTransfer.strains = BMat*solution + BMatWarp2*warpSol;
  //  interfaceElement->getMaterialFormulation(pointers, i)->getMaterialData(pointers, materialTransfer,i);
  //  auto coor = interfaceElement->getLocalCoordinate(pointers,i);
  //  std::cout << i.xi << ";"<<i.eta<<";" << i.zeta << ";" << coor.transpose().format(CleanFmt) << " " << materialTransfer.stresses.transpose().format(CleanFmt) << std::endl;
  //}

}

auto EL302_BeamCoupling3D::getBMatrixLinear(PointerCollection& pointers, Geometry::H1Shapes &shapes, IntegrationPoint &intPoint, FiniteElement::beamInterfaceElement3D &elem) -> Types::Matrix6X<prec>
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
  
  for (auto i=0;i<numShapes-2;++i){
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

  indexType pos = numShapes-2;
  pos *= 3;
  BMat.block(0, pos, 1, 3) = shapes.shapeDeriv(0, numShapes - 2) * A1;
  BMat.block(3, pos, 1, 3) = shapes.shapeDeriv(0, numShapes - 2) * A2;
  BMat.block(4, pos, 1, 3) = shapes.shapeDeriv(0, numShapes - 2) * A3;

  pos += 3;
  BMat.block(0, pos, 1, 3) = shapes.shapeDeriv(0,numShapes-1)*(localCoor(2)*A2-localCoor(1)*A3);

  BMat.block(3, pos, 1, 3) = -shapes.shapeDeriv(0,numShapes-1)*localCoor(2)*A1;
  BMat.block(3, pos, 1, 3) += -shapes.shapes(numShapes-1)*A3;

  BMat.block(4, pos, 1, 3) = shapes.shapeDeriv(0,numShapes-1)*localCoor(1)*A1;
  BMat.block(4, pos, 1, 3) += shapes.shapes(numShapes-1)*A2;
  

  return BMat;
}



} // namespace HierAMuS::Elementformulations
