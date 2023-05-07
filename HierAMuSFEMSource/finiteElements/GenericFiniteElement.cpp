// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include <iostream>
#include <pointercollection/pointercollection.h>


#include <finiteElements/GenericFiniteElement.h>

#include <elementFormulations/GenericElementFormulation.h>

#include <geometry/GeometryData.h>

#include <materials/Material.h>

#include <finiteElements/ElementTypes.h>

#include <solver/GenericSolutionState.h>

#include <Eigen/Dense>
#include <sstream>

#include "materials/GenericMaterialFormulation.h"

namespace HierAMuS::FiniteElement {

GenericFiniteElement::GenericFiniteElement() { this->material = nullptr; }

GenericFiniteElement::~GenericFiniteElement() = default;



void GenericFiniteElement::setAllNodeBoundaryConditionMeshId(ptrCol &pointers,
                                                             indexType meshId,
                                                             indexType dof) {}

void GenericFiniteElement::insertStiffnessResidual(
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, 1, Eigen::Dynamic> &residual,
    std::vector<indexType> &eqIds, std::vector<dofStatus> &eqStatus) {}

void GenericFiniteElement::GenericSetDegreesOfFreedom(PointerCollection& pointers) {

  this->material->getElementFormulation(pointers)->setDegreesOfFreedom(pointers,this);
}

void GenericFiniteElement::GenericAdditionalOperations(PointerCollection& pointers) {
  this->material->getElementFormulation(pointers)->AdditionalOperations(pointers, this);
}

auto GenericFiniteElement::getDofs(PointerCollection& pointers) -> std::vector<DegreeOfFreedom *> {
  return this->material->getElementFormulation(pointers)->getDofs(pointers, this);
}

void GenericFiniteElement::GenericSetTangentResidual(
  PointerCollection& pointers,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  this->material->getElementFormulation(pointers)->setTangentResidual(pointers, this,
                                                              stiffness, residual, Dofs);
}

auto GenericFiniteElement::computeNorm(GenericFiniteElement::ptrCol& pointers, NormTypes type) -> prec {
  return this->material->getElementFormulation(pointers)->computeNorm(pointers, this, type);
}


void GenericFiniteElement::GenericSetMass(
  PointerCollection& pointers,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {
  this->material->getElementFormulation(pointers)->setMass(pointers, this, stiffness, residual,
                                                   Dofs);
}

void GenericFiniteElement::getElementsLocalNodalReactions(
    PointerCollection &ptrCol, std::map<indexType, std::vector<prec>> &vReacs) {
  this->material->getElementFormulation(ptrCol)->getElementsLocalNodalReactions(
      ptrCol, this, vReacs);
}

void GenericFiniteElement::getSolution(ptrCol &pointers,
                                       std::vector<DegreeOfFreedom *> &Dofs,
                                       Types::VectorX<prec> &solution) {

  auto solstate = pointers.getSolutionState();

  solution = solstate->getSolution(Dofs);
}

auto GenericFiniteElement::getSolution(ptrCol &pointers,
                                       std::vector<DegreeOfFreedom *> &Dofs) -> Types::VectorX<prec> {

  auto solstate = pointers.getSolutionState();

  return solstate->getSolution(Dofs);
}

void GenericFiniteElement::getVelocity(ptrCol &pointers,
                                       std::vector<DegreeOfFreedom *> &Dofs,
                                       Types::VectorX<prec> &solution) {

  auto solstate = pointers.getSolutionState();

  solution = solstate->getVelocity(Dofs);
}

void GenericFiniteElement::getAcceleration(ptrCol &pointers,
                                           std::vector<DegreeOfFreedom *> &Dofs,
                                           Types::VectorX<prec> &solution) {

  auto solstate = pointers.getSolutionState();

  solution = solstate->getAcceleration(Dofs);
}

void GenericFiniteElement::setParaviewCellData(ptrCol &pointers,
                                               vtkPlotInterface &catalyst) {}

auto GenericFiniteElement::getNumberOfIntegrationPoints(PointerCollection& pointers) -> indexType {
  return this->material->getElementFormulation(pointers)->getNumberOfIntergrationPoints(pointers, this);
}

auto GenericFiniteElement::getElementHistoryDataStructure(PointerCollection& pointers)
-> const HistoryDataStructure & {
  return this->material->getElementFormulation(pointers)->getHistoryDataStructure();
}

auto GenericFiniteElement::getMaterialHistoryDataStructure(PointerCollection& pointers)
-> const HistoryDataStructure & {
  return this->material->getMaterialFormulation(pointers)->getHistoryDataStructure(pointers);
}

auto GenericFiniteElement::getHistoryDataSetUp(PointerCollection& pointers) -> HistoryDataSetup
{
  auto elemHistStruct = this->getElementHistoryDataStructure(pointers);
  auto matHistStrcut = this->getMaterialHistoryDataStructure(pointers);
  indexType numGp =
      this->getNumberOfIntegrationPoints(pointers);
  HistoryDataSetup setup;
  setup.setNumberOfIntegrationPoints(numGp);
  setup.setElementId(this->id);
  setup.setNHistPerIntPointElement(elemHistStruct.getNumberOfUpdateValues(),
                                   elemHistStruct.getNumberOfConstValues());
  setup.setNHistPerIntPointMaterial(matHistStrcut.getNumberOfUpdateValues(),
                                    matHistStrcut.getNumberOfConstValues());
  return setup;
}

auto GenericFiniteElement::getHistoryDataIterator(ptrCol &pointers)
    -> HistoryDataIterator {
  auto histIt = pointers.getSolutionState()->getHistoryData(this->id);
  histIt.setElementStructure(this->getElementHistoryDataStructure(pointers));
  histIt.setMaterialStructure(this->getMaterialHistoryDataStructure(pointers));
  return histIt;
}

void GenericFiniteElement::updateRVEHistory(PointerCollection &pointers)
{
  this->material->getElementFormulation(pointers)->updateRVEHistory(pointers,
                                                                       this);
}

void GenericFiniteElement::toParaviewAdapter(
    GenericFiniteElement::ptrCol &pointers, vtkPlotInterface &catalyst,
    ParaviewSwitch ToDo) {
  this->material->getElementFormulation(pointers)->toParaviewAdaper(pointers, this,
                                                            catalyst, ToDo);
}

auto GenericFiniteElement::getVertex(GenericFiniteElement::ptrCol &pointers,
                                     indexType localNumber)
    -> Geometry::Vertex & {

  throw std::runtime_error("Method getVertex not implemented to the element");
  return pointers.getGeometryData()->getVertex(1);
}

auto GenericFiniteElement::getEdge(ptrCol &pointers, indexType localNumber)
    -> Geometry::Edges & {
  throw std::runtime_error("Method getEdge not implemented for Element");
  return pointers.getGeometryData()->getEdge(1);
};

auto GenericFiniteElement::getIntegrationPoints(ptrCol &pointers)
    -> IntegrationPoints {
  auto temp = HierAMuS::FiniteElement::GenericFiniteElement::ptrCol::getIntegrationPoints(this->id);
  return temp;
}

} /* namespace HierAMuS */
