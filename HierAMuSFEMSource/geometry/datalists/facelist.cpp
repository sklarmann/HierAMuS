// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "geometry/GeometryTypes.h"
#include <geometry/datalists/facelist.h>

#include <sstream>

namespace HierAMuS::Geometry {

facelist::facelist() {
  this->maxFace = 0;
  this->numFaces = 0;
}

facelist::~facelist() {
  this->linearTriangles.clear();
  this->linearQuadrilaterals.clear();
  this->quadraticQuadrilaterals.clear();
  this->faceIdMap.clear();
  this->maxFace = 0;
  this->numFaces = 0;
}

auto facelist::newElement(GeometryTypes type, indexType number) -> indexType {

  if (this->faceIdMap.find(number) == this->faceIdMap.end()) {
    switch (type) {
    case GeometryTypes::LinearTriangle: {
      indexType pos1 = this->linearTriangles.size();
      this->faceIdMap[number] = {type, pos1};
      this->linearTriangles.emplace_back();
      this->linearTriangles[pos1].setId(number);
    } break;
    case GeometryTypes::LinearQuadrilateral: {
      indexType pos1 = this->linearQuadrilaterals.size();
      this->faceIdMap[number] = {type, pos1};
      this->linearQuadrilaterals.emplace_back();
      this->linearQuadrilaterals[pos1].setId(number);

    } break;
    case GeometryTypes::QuadraticQuadrilateral: {
      indexType pos1 = this->quadraticQuadrilaterals.size();
      this->faceIdMap[number] = {type, pos1};
      this->quadraticQuadrilaterals.emplace_back();
      this->quadraticQuadrilaterals[pos1].setId(number);

    } break;
    case GeometryTypes::ScaledBoundary2D: {
      indexType pos1 = this->scaledBoundary.size();
      this->faceIdMap[number] = {type, pos1};
      this->scaledBoundary.emplace_back();
      this->scaledBoundary[pos1].setId(number);
    } break;
    default:
      std::stringstream msg;
      msg << "Trying to add an face with the wrong type!";
      throw std::runtime_error(msg.str());
    }
    if (number > this->maxFace) {
      this->maxFace = number;
    }
    ++this->numFaces;
  }
  return number;
}

auto facelist::getFace(indexType number) -> Faces & {
  if (this->faceIdMap.find(number) == this->faceIdMap.end()) {
    std::stringstream msg;
    msg << "Requested Edge with number " << number << " does not exist!";
    throw std::runtime_error(msg.str());
  }
  auto TypePos = this->faceIdMap[number];

  switch (TypePos.first) {
  case GeometryTypes::LinearTriangle:
    return this->linearTriangles[TypePos.second];
  case GeometryTypes::LinearQuadrilateral:
    return this->linearQuadrilaterals[TypePos.second];
  case GeometryTypes::QuadraticQuadrilateral:
    return this->quadraticQuadrilaterals[TypePos.second];
  case GeometryTypes::ScaledBoundary2D:
    return this->scaledBoundary[TypePos.second];
  default:
    std::stringstream msg;
    msg << "Requested Face with number " << number << " does not exist!"
        << std::endl;
    throw std::runtime_error(msg.str());
  }
}

auto facelist::getItemLocalNumber(indexType number) -> Faces * {
  indexType posA = this->linearTriangles.size();
  indexType posB = posA + this->linearQuadrilaterals.size();
  indexType posC = posB + this->scaledBoundary.size();
  indexType posD = posC + this->quadraticQuadrilaterals.size();
  if (number < posA) {
    return &this->linearTriangles[number];
  }
  else if (number < posB) {
    return &this->linearQuadrilaterals[number - posA];
  }
  else if (number < posC) {
    return &this->scaledBoundary[number - posB];
  }
  else if (number < posD) {
    return &this->quadraticQuadrilaterals[number - posC];
  }
  return nullptr;
}

auto facelist::getFaceNumberByVertexNumbers(indexType vert1, indexType vert2, indexType vert3) -> indexType{

  for(auto &face : this->linearTriangles){
    if(face.hasVertices(vert1, vert2, vert3)){
      return face.getId();
    }
  }
  for(auto &face : this->linearQuadrilaterals){
    if(face.hasVertices(vert1, vert2, vert3)){
      return face.getId();
    }
  }
  for(auto &face : this->scaledBoundary){
    if(face.hasVertices(vert1, vert2, vert3)){
      return face.getId();
    }
  }

  return -1;
}


} // namespace HierAMuS::Geometry
