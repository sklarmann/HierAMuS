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

GenericFiniteElement::GenericFiniteElement() { this->m_material = nullptr; }

GenericFiniteElement::~GenericFiniteElement() = default;



void GenericFiniteElement::setAllNodeBoundaryConditionMeshId(ptrCol &pointers,
                                                             indexType meshId,
                                                             indexType dof) {}
 



auto GenericFiniteElement::getDofs(PointerCollection& pointers) -> std::vector<DegreeOfFreedom *> {
  return this->m_material->getElementFormulation(pointers)->getDofs(pointers, this);
}



auto GenericFiniteElement::computeNorm(GenericFiniteElement::ptrCol& pointers, NormTypes type) -> prec {
  return this->m_material->getElementFormulation(pointers)->computeNorm(pointers, this, type);
}


void GenericFiniteElement::GenericSetMass(
  PointerCollection& pointers,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {
  this->m_material->getElementFormulation(pointers)->setMass(pointers, this, stiffness, residual,
                                                   Dofs);
}

void GenericFiniteElement::getElementsLocalNodalReactions(
    PointerCollection &ptrCol, std::map<indexType, std::vector<prec>> &vReacs) {
  this->m_material->getElementFormulation(ptrCol)->getElementsLocalNodalReactions(
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

  return pointers.getSolutionState()->getSolution(Dofs);
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


auto GenericFiniteElement::getNumberOfIntegrationPoints(PointerCollection& pointers) -> indexType {
  return this->m_material->getElementFormulation(pointers)->getNumberOfIntergrationPoints(pointers, this);
}

auto GenericFiniteElement::getElementHistoryDataStructure(PointerCollection& pointers)
-> const HistoryDataStructure & {
  return this->m_material->getElementFormulation(pointers)->getHistoryDataStructure();
}

auto GenericFiniteElement::getMaterialHistoryDataStructure(PointerCollection& pointers)
-> const HistoryDataStructure & {
  return this->m_material->getMaterialFormulation(pointers)->getHistoryDataStructure(pointers);
}

auto GenericFiniteElement::getHistoryDataSetUp(PointerCollection& pointers) -> HistoryDataSetup
{
  auto elemHistStruct = this->getElementHistoryDataStructure(pointers);
  auto matHistStrcut = this->getMaterialHistoryDataStructure(pointers);
  indexType numGp =
      this->getNumberOfIntegrationPoints(pointers);
  HistoryDataSetup setup;
  setup.setNumberOfIntegrationPoints(numGp);
  setup.setElementId(this->m_id);
  setup.setNHistPerIntPointElement(elemHistStruct.getNumberOfUpdateValues(),
                                   elemHistStruct.getNumberOfConstValues());
  setup.setNHistPerIntPointMaterial(matHistStrcut.getNumberOfUpdateValues(),
                                    matHistStrcut.getNumberOfConstValues());
  return setup;
}

auto GenericFiniteElement::getHistoryDataIterator(ptrCol &pointers)
    -> HistoryDataIterator {
  auto histIt = pointers.getSolutionState()->getHistoryData(this->m_id);
  histIt.setElementStructure(this->getElementHistoryDataStructure(pointers));
  histIt.setMaterialStructure(this->getMaterialHistoryDataStructure(pointers));
  return histIt;
}

void GenericFiniteElement::updateRVEHistory(PointerCollection &pointers)
{
  this->m_material->getElementFormulation(pointers)->updateRVEHistory(pointers,
                                                                       this);
}


auto GenericFiniteElement::getVertex(GenericFiniteElement::ptrCol &pointers,
                                     indexType localNumber)
    -> Geometry::VertexData & {

  throw std::runtime_error("Method getVertex not implemented to the element");
  return pointers.getGeometryData()->getVertexData(1);
}

auto GenericFiniteElement::getEdge(ptrCol &pointers, indexType localNumber)
    -> Geometry::EdgesData & {
  throw std::runtime_error("Method getEdge not implemented for Element");
  return pointers.getGeometryData()->getEdgeData(1);
};

auto GenericFiniteElement::getIntegrationPoints(ptrCol &pointers)
    -> IntegrationPoints {
  return IntegrationPointsManagement::getIntegrationsPoints(this->m_id);
}

auto GenericFiniteElement::getMaterialFormulation(PointerCollection &pointers)
    -> std::shared_ptr<HierAMuS::Materials::GenericMaterialFormulation> {
  return this->m_material->getMaterialFormulation(pointers);
}
auto GenericFiniteElement::getMaterialFormulation(PointerCollection &pointers,
                                                  IntegrationPoint &ip)
    -> std::shared_ptr<HierAMuS::Materials::GenericMaterialFormulation> {
  return this->m_material->getMaterialFormulation(pointers);
};
auto GenericFiniteElement::getMaterialId() -> indexType {
  return this->m_material->getNumber();
}

void GenericFiniteElement::request_element_data_field(
    PointerCollection &pointers, indexType fieldId, indexType rows,
    indexType cols) {
  pointers.getSolutionState()->request_element_data_field(this->m_id, fieldId,
                                                          rows, cols);
}

auto GenericFiniteElement::get_element_data_field(PointerCollection &pointers,
                                                  indexType fieldId)
    -> Types::MatrixXX<prec> & {
  return pointers.getSolutionState()->get_element_data_field(m_id, fieldId);
}

void GenericFiniteElement::set_element_data_field(PointerCollection &pointers,
                                                  indexType fieldId,
                                                  Types::MatrixXX<prec> &data) {
  pointers.getSolutionState()->set_element_data_field(m_id, fieldId, data);
}

auto GenericFiniteElement::getIncrementalSolution(
    ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs)
    -> Types::VectorX<prec> {
  return pointers.getSolutionState()->getIncrementalSolution(Dofs);
}

auto GenericFiniteElement::getNewtonSolution(
    ptrCol &pointers, std::vector<DegreeOfFreedom *> &Dofs)
    -> Types::VectorX<prec> {
  return pointers.getSolutionState()->getNewtonSolution(Dofs);
}

} /* namespace HierAMuS */
