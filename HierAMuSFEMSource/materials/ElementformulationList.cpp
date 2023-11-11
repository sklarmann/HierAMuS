// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <materials/ElementformulationList.h>
#include <pointercollection/pointercollection.h>


#include <elementFormulations/EL101_Bernoulli2D.h>
#include <elementFormulations/EL102_Timoshenko2D.h>
#include <elementFormulations/EL103_Timoshenko3D.h>
#include <elementFormulations/EL104_TimoshenkoPrism.h>
#include <elementFormulations/LSFEMBernoulli.h>

#include <elementFormulations/EL201_2DShell.h>
#include <elementFormulations/EL202_Piezo2D.h>
#include "elementFormulations/EL290_PythonElement.h"

#include <elementFormulations/EL203_BeamInterface2D.h>
#include <elementFormulations/EL204_BeamInterface2D.h>
#include <elementFormulations/EL205_HDivTest.h>
#include <elementFormulations/EL206_Plate.h>
#include <elementFormulations/EL207_FaceConstraint.h>

#include <elementFormulations/EL300_Solid3DLinear.h>
#include <elementFormulations/EL301_Piezo3DLinear.h>
#include "elementFormulations/EL302_BeamCoupling3D.h"
#include "elementFormulations/EL303_ThermoMechanikSolid3D.h"
#include "elementFormulations/EL304_QPVolumeElement.h"
#include "elementFormulations/EL307_VolumeConstraint.h"

#include "control/ParameterList.h"


namespace HierAMuS::Materials {

ElementFormulationList::ElementFormulationList() = default;

ElementFormulationList::~ElementFormulationList() = default;

void ElementFormulationList::addElementFormulation(
    PointerCollection &pointers, indexType number, indexType elementFormulation,
    ParameterList &elementparameters) {
  while (number + 1 > static_cast<indexType>(this->Elements.size())) {
    this->Elements.push_back(nullptr);
  }

  switch (elementFormulation) {

    // Linear 2D Bernoulli Beam
  case 101:
    this->Elements[number] =
        std::make_shared<Elementformulations::EL101_Bernoulli2D>(&pointers);
    break;
    // Finite Rotation 2D Timoshenko Beam
  case 102:
    this->Elements[number] =
        std::make_shared<Elementformulations::EL102_Timoshenko2D>(&pointers);
    break;
    // Least squares Bernoulli Beam
  case 103:
    this->Elements[number] =
        std::make_shared<Elementformulations::EL103_Timoshenko3D>(&pointers);
    break;
    // 3D Timoshenko Beam with Prisms
  case 104:
    this->Elements[number] =
        std::make_shared<Elementformulations::EL104_TimoshenkoPrism>(&pointers);
    break;
    // Least squares Bernoulli Beam
  case 125:
    this->Elements[number] =
        std::make_shared<Elementformulations::LSFEMBernoulli>(&pointers);
    break;
    // 2D Shell Element
  case 201:
    this->Elements[number] =
        std::make_shared<Elementformulations::EL201_2DShell>(&pointers);
    break;
    // 2D Piezo Element
  case 202:
    this->Elements[number] =
        std::make_shared<Elementformulations::EL202_Piezo2D>(&pointers);
    break;
    // 2D Beam Interface Element
  case 203:
    this->Elements[number] =
        std::make_shared<Elementformulations::EL203_BeamInterface2D>(&pointers);
    break;
  case 204:
    this->Elements[number] =
        std::make_shared<Elementformulations::EL204_BeamInterface2D>(&pointers);
    break;
  case 205:
    this->Elements[number] =
		std::make_shared<Elementformulations::EL205_HDivTest>(&pointers);
    break;
  case 206:
    this->Elements[number] =
		std::make_shared<Elementformulations::EL206_Plate>(&pointers);
    break;
  case 207:
    this->Elements[number] =
		std::make_shared<Elementformulations::EL207_FaceConstraint>(&pointers);
    break;
  case 290:
    this->Elements[number] =
        std::make_shared<Elementformulations::EL290_2DPythonElement>(&pointers);
    // 3D Brick Element
  case 300:
    this->Elements[number] =
        std::make_shared<Elementformulations::EL300_Solid3DLinear>(&pointers);
    break;
    // 3D Piezo Brick Element
  case 301:
    this->Elements[number] =
        std::make_shared<Elementformulations::EL301_Piezo3DLinear>(&pointers);
    break;
  case 302:
    this->Elements[number] =
        std::make_shared<Elementformulations::EL302_BeamCoupling3D>(&pointers);
    break;
  case 303:
    this->Elements[number] =
        std::make_shared<Elementformulations::EL303_ThermoMechanikSolid3D>(
            &pointers);
    break;
  case 304:
    this->Elements[number] =
        std::make_shared<Elementformulations::EL304_QPVolumeElement>(
            &pointers);
    break;
  case 307:
    this->Elements[number] =
        std::make_shared<Elementformulations::EL307_VolumeConstraint>(
            &pointers);
    break;
  default:
    std::string temp;
    temp = "Selected Elementformulation does not exist!";
    throw std::runtime_error(temp);
  }
  auto temp = this->Elements[number];
  temp->readData(pointers,elementparameters);
}

std::shared_ptr<Elementformulations::GenericElementFormulation>
ElementFormulationList::getElementFormulation(indexType number) {
  return this->Elements[number];
}
void ElementFormulationList::addElementFormulation(indexType number,
    std::shared_ptr<Elementformulations::GenericElementFormulation> element) {
  while (number + 1 > static_cast<indexType>(this->Elements.size())) {
    this->Elements.push_back(nullptr);
  }
  this->Elements[number] = element;
}
} // namespace HierAMuS
