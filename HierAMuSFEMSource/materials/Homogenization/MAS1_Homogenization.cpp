// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include "MAS1_Homogenization.h"

#include "solver/StaticSolutionStateHomogenization.h"

#include "control/ParameterList.h"


#include <string>

#include "control/HandlingStructs.h"
#include "pointercollection/pointercollection.h"

namespace HierAMuS::Materials {

MAS1_Homogenization::MAS1_Homogenization(PointerCollection *ptrCol)
    : GenericMaterialFormulation(ptrCol), m_maxIterations(10), m_refNorm(0), m_refENorm(0)
{
#ifdef OPENMP
  omp_init_lock(&m_normLock);
#endif
}

MAS1_Homogenization::~MAS1_Homogenization() {
#ifdef OPENMP
  omp_destroy_lock(&m_normLock);
#endif
};

void MAS1_Homogenization::readData(PointerCollection& pointers, ParameterList &list) {
  pointers.setRVEs();
  m_maxIterations = list.getIndexVal("maxIterations");
}

void MAS1_Homogenization::getMaterialData(
  PointerCollection& pointers, MaterialTransferData &materialInOut, IntegrationPoint& ip) {

  // get Copy of the m_RVE
  auto RVEcp = this->m_RVE->getShallowCopy();
  auto fname = this->getFileName(RVEcp, ip);

  // read currents m_RVE data
  RVEcp->RVEDataFromFile(fname);
  

  // set current strains to m_RVE
  auto solstate = std::dynamic_pointer_cast<StaticSolutionStateHomogenization>(
      RVEcp->getSolutionState());
  solstate->setStrains(*RVEcp, materialInOut.strains);
  

  // perform first solution of the system of equations
  solstate->assembleSystem(*RVEcp);
  solstate->factorize();
  solstate->solve(*RVEcp);
  solstate->updateSolution(*RVEcp);

  prec en0 = solstate->energyNorm();
  prec res0 = solstate->residual();
  prec strainNorm = materialInOut.strains.norm();

  indexType iteration = 1;
  prec enorm = prec(1);
  prec resnorm = prec(1);
  if (strainNorm == prec(0))
  {
    iteration = m_maxIterations;
  } else
  {
#ifdef OPENMP
    omp_set_lock(&m_normLock);
#endif

    if (abs(en0) > m_refENorm)
    {
      m_refENorm = abs(en0);
    }
    if (abs(res0) > m_refNorm)
    {
      m_refNorm = abs(res0);
    }
    enorm = abs(en0 / m_refENorm);
    resnorm = abs(res0 / m_refNorm);
#ifdef OPENMP
    omp_unset_lock(&m_normLock);
#endif
  }

  while (iteration < m_maxIterations &&
         abs(enorm) > std::numeric_limits<prec>::epsilon() * prec(1000) &&
         abs(resnorm) > std::numeric_limits<prec>::epsilon() * prec(1000)) {
    solstate->assembleSystem(*RVEcp);
    solstate->factorize();
    solstate->solve(*RVEcp);
    solstate->updateSolution(*RVEcp);

    enorm = solstate->energyNorm() / m_refENorm;
    resnorm = solstate->residual() / m_refNorm;
    ++iteration;
  }
  
  solstate->homogenize(*RVEcp);

  auto homData = solstate->getHomogenizedData();
  materialInOut.materialTangent = homData.C;
  materialInOut.stresses = homData.sigma;

  auto &RVELogger = RVEcp->getSPDLogger();
  RVELogger.info("\nRVE converged after: {:>4} / {:<4} iterations."
                "With energy norm: {:12.6e}\n"
                "Residual norm:   {:12.6e}\n",
                iteration, m_maxIterations, enorm, resnorm);


  
  RVEcp->RVEDataToFile(fname);
}

void MAS1_Homogenization::initRVE(PointerCollection &pointers,
                                  IntegrationPoint &ip)
{
  auto RVEcp = this->m_RVE->getShallowCopy();
  auto fname = this->getFileName(RVEcp, ip);
  RVEcp->RVEDataToFile(fname);
}

void MAS1_Homogenization::setRVE(PointerCollection &pointers,
                                 PointerCollection &RVE) {
  this->m_RVE = &RVE;
}

void MAS1_Homogenization::updateRVEHistory(PointerCollection &pointers,
                                           IntegrationPoint &ip)
{
  auto RVEcp = this->m_RVE->getShallowCopy();
  auto fname = this->getFileName(RVEcp, ip);

  RVEcp->RVEDataFromFile(fname);

  RVEcp->getSolutionState()->nextSolutionStep();

  RVEcp->RVEDataToFile(fname);
}

auto MAS1_Homogenization::getFileName(
    std::shared_ptr<PointerCollection> pointers,
                                      IntegrationPoint &ip)
    -> std::string {
  std::stringstream fname;
  fname << pointers->getInfoData()->fileNames[FileHandling::infile] << "."
        << ip.elementId << "." << ip.sectionNumber << "." << ip.gpNumber
        << ".sol";
  return fname.str();
}

} // namespace HierAMuS
