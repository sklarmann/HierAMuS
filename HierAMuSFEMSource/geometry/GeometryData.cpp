// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <cstddef>

#include <geometry/GeometryData.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <pointercollection/pointercollection.h>

#include <geometry/BeamInterface2D.h>
#include <geometry/BeamInterface3D.h>
#include <geometry/Edges.h>
#include <geometry/Faces.h>
#include <geometry/LinearBrick.h>
#include <geometry/LinearEdge.h>
#include <geometry/LinearPrism.h>
#include <geometry/LinearQuadrilateral.h>
#include <geometry/LinearTriangle.h>
#include <geometry/Special.h>
#include <geometry/Vertex.h>
#include <geometry/Volumes.h>

#include <math/Plane.h>

#include <algorithm>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <vector>

namespace HierAMuS::Geometry {

GeometryData::GeometryData() {}

GeometryData::~GeometryData() {}

auto GeometryData::requestNewVert() -> indexType {
  indexType number = this->vertices.lastElement() + 1;
  // indexType number = this->vertList.getLastElementsIndex() + 1;
  this->vertices.newElement(GeometryTypes::Vertex, number);
  this->vertices.getVertex(number).setId(number);
  return number;
}

auto GeometryData::getVertex(indexType num) -> Vertex & {
  return this->vertices.getVertex(num);
}

auto GeometryData::getEdge(indexType num) -> Edges & {
  return this->edges.getEdge(num);
}

auto GeometryData::getFace(indexType num) -> Faces * {
  return &this->faces.getFace(num);
}

auto GeometryData::getVolume(indexType num) -> Volumes * {
  return &this->volumes.getVolume(num);
}

auto GeometryData::getSpecial(indexType num) -> Special * {
  Special *ret = nullptr;
  if (this->specialList.inList(num)) {
    ret = this->specialList[num];
  } else {
    std::stringstream msg;
    msg << "In GeometryData->getSpecial(), requested Special element " << num
        << " is not in list!";
    throw std::runtime_error(msg.str());
  }
  return ret;
}

auto GeometryData::getVertexClosestTo(Types::Vector3<prec> point) -> Vertex & {
  return this->vertices.getVertexClosestTo(point);
}

auto GeometryData::requestNewGeometryObject(GeometryTypes type) -> indexType {
  indexType ret = -1;
  switch (type) {
  case GeometryTypes::Vertex:
    ret = this->vertices.lastElement();
    break;
  case GeometryTypes::LinearEdge:
  case GeometryTypes::QuadraticEdge:
    ret = this->edges.lastElement();
    break;
  case GeometryTypes::LinearTriangle:
  case GeometryTypes::LinearQuadrilateral:
  case GeometryTypes::QuadraticQuadrilateral:
  case GeometryTypes::ScaledBoundary2D:
    ret = this->faces.lastElement();
    break;
  case GeometryTypes::LinearBrick:
  case GeometryTypes::LinearPrism:
    ret = this->volumes.lastElement();
    break;
  case GeometryTypes::BeamInterface2D:
  case GeometryTypes::BeamInterface3D:
    ret = this->specialList.getLastElementsIndex();
    break;
  default:
    throw std::invalid_argument("Selected Geometry type to add does not exist");
    break;
  }
  ++ret;
  this->requestNewGeometryObject(type, ret);
  return ret;
}

auto GeometryData::requestNewGeometryObject(GeometryTypes type,
                                            indexType number) -> indexType {
  indexType ret = 0;
  switch (type) {
  case GeometryTypes::Vertex:
    this->vertices.newElement(type, number);
    this->vertices.getVertex(number).setId(number);
    // this->vertList.insertAtIndex(number, new Vertex);
    // this->vertList[number]->setId(number);
    ret = number;
    break;
  case GeometryTypes::LinearEdge:
  case GeometryTypes::QuadraticEdge:
    this->edges.newElement(type, number);
    ret = number;
    ;
    break;
  case GeometryTypes::LinearTriangle:
  case GeometryTypes::LinearQuadrilateral:
  case GeometryTypes::QuadraticQuadrilateral:
  case GeometryTypes::ScaledBoundary2D:
    this->faces.newElement(type, number);
    ret = number;
    break;
  case GeometryTypes::LinearBrick:
  case GeometryTypes::LinearPrism:
    this->volumes.newElement(type, number);
    ret = number;
    break;
  case GeometryTypes::BeamInterface2D:
    this->specialList.insertAtIndex(number, new BeamInterface2D);
    this->specialList[number]->setId(number);
    ret = number;
    break;
  case GeometryTypes::BeamInterface3D:
    this->specialList.insertAtIndex(number, new BeamInterface3D);
    this->specialList[number]->setId(number);
    ret = number;
    break;
  default:
    // TODO throw exception
    throw std::runtime_error("In GeometryData->requstNewGeometryObject(), "
                             "tried to add unknown geometry type!");
    break;
  }
  return ret;
}

auto GeometryData::getGeometryElement(GeometryTypes type, indexType num)
    -> Base * {

  switch (type) {
  case GeometryTypes::Vertex:
    return &this->vertices.getVertex(num);
    break;
  case GeometryTypes::Edges:
    return &this->edges.getEdge(num);
    break;
  case GeometryTypes::Faces:
    return &this->faces.getFace(num);
    break;
  case GeometryTypes::Volumes:
    return &this->volumes.getVolume(num);
    break;
  case GeometryTypes::Special:
    if (!this->specialList.inList(num)) {
      std::stringstream msg;
      msg << "Error in GeometryData->getGeometryElement(): Requested Volume "
          << num << " not in list!";
      throw std::runtime_error(msg.str());
    }
    return this->specialList[num];
    break;
  default:
    // TODO throw exception
    std::cout << "getGeometryElement Error" << std::endl;
  }
  return nullptr;
}

void GeometryData::getGeometricElementInPlane(std::vector<prec> normal,
                                              std::vector<prec> point,
                                              GeometryTypes type,
                                              std::stack<Base *> &elems) {

  Plane plane;
  plane.set(normal, point);
  std::vector<indexType> vertNums;

  for (auto vertice : this->vertices) {
    if (plane.inPlane(vertice->getCoordinates())) {
      vertNums.push_back(vertice->getId());
    }
  }

  switch (type) {
  case GeometryTypes::Vertex: {
    for (auto &vertNum : vertNums) {
      // elems.push(this->vertList[*it]);
      elems.push(&this->vertices.getVertex(vertNum));
    }
  } break;
  case GeometryTypes::Edges: {
    if (!vertNums.empty()) {
      std::vector<indexType> edgeVerts;
      for (auto edge : this->edges) {
        bool add = true;
        edgeVerts.clear();
        edge->getVerts(edgeVerts);
        for (auto &edgeVert : edgeVerts) {
          if (std::find(vertNums.begin(), vertNums.end(), edgeVert) ==
              vertNums.end()) {
            add = false;
          }
        }
        if (add) {
          elems.push(edge);
        }
      }
    }
    break;
  case GeometryTypes::Faces: {
    if (!vertNums.empty()) {
      std::vector<indexType> faceVerts;
      for (auto face : this->faces) {
        faceVerts.clear();
        face->getVerts(faceVerts);
        bool add = true;
        for (auto &faceVert : faceVerts) {
          if (std::find(vertNums.begin(), vertNums.end(), faceVert) ==
              vertNums.end()) {
            add = false;
          }
        }
        if (add) {
          elems.push(face);
        }
      }
    }
  }
  } break;
  default:
    break;
  }
}

auto GeometryData::getVerticesInPlane(PointerCollection &pointers,
                                      const Types::Vector3<prec> &normal,
                                      const Types::Vector3<prec> &point)
    -> std::vector<Geometry::Vertex *> {
  Plane plane;
  plane.set(normal, point);

  std::vector<Vertex *> verts;
  for (auto v : this->vertices) {
    if (plane.inPlane(v->getCoordinates()))
      verts.push_back(v);
  }
  return verts;
}

auto GeometryData::getEdgesInPlane(PointerCollection &pointers,
                                   const Types::Vector3<prec> &normal,
                                   const Types::Vector3<prec> &point)
    -> std::vector<Geometry::Edges *> {
  Plane plane;
  plane.set(normal, point);

  auto verts = this->getVerticesInPlane(pointers, normal, point);

  std::map<indexType, Edges *> idEdgeMap;
  IntegrationPoint ip;
  ip.xi = 0;
  for (auto v : verts) {
    auto pedges = v->getConnectedEdges();
    for (auto edgenum : pedges) {
      auto &edge = this->getEdge(edgenum);
      if (plane.inPlane(edge.getCoordinates(pointers, ip)))
        idEdgeMap[edge.getId()] = &edge;
    }
  }

  std::vector<Geometry::Edges *> edgeList;
  for (auto j : idEdgeMap) {
    edgeList.push_back(j.second);
  }

  return edgeList;
}

auto GeometryData::getFacesInPlane(PointerCollection &pointers,
                                   const Types::Vector3<prec> &normal,
                                   const Types::Vector3<prec> &point)
    -> std::vector<Geometry::Faces *> {
  Plane plane;
  plane.set(normal, point);

  auto verts = this->getVerticesInPlane(pointers, normal, point);

  std::map<indexType, Faces *> idFaceMap;
  IntegrationPoint ip;
  ip.xi = 0;
  ip.eta = 0;
  for (auto v : verts) {
    auto pfaces = v->getConnectedFaces();
    for (auto facenum : pfaces) {
      auto face = this->getFace(facenum);
      if (plane.inPlane(face->getCoordinates(pointers, ip)))
        idFaceMap[face->getId()] = face;
    }
  }

  std::vector<Geometry::Faces *> faceList;
  for (auto j : idFaceMap) {
    faceList.push_back(j.second);
  }

  return faceList;
}

void GeometryData::print(PointerCollection &pointers) {

  auto &Logger = pointers.getSPDLogger();

  Logger.info(
      "\nGeometry Data informations:\n"
      "   Number of vertices:   {:>12}\n"
      "   Number of edges:      {:>12}\n"
      "   Number of faces:      {:>12}\n"
      "   Number of volumes:    {:>12}",
      this->vertices.getNumberOfElements(), this->edges.getNumberOfElements(),
      this->faces.getNumberOfElements(), this->volumes.getNumberOfElements());

  for (auto vertice : this->vertices) {
    vertice->print(pointers);
    Logger.debug("{:-<100}", "");
  }
  for (auto edge : this->edges) {
    edge->print(pointers);
    Logger.debug("{:-<100}", "");
  }
  for (auto face : this->faces) {
    face->print(pointers);
    Logger.debug("{:-<100}", "");
  }
  for (auto vol : this->volumes) {
    vol->print(pointers);
    Logger.debug("{:-<100}", "");
  }
}

auto GeometryData::getEdgeNumberByVerts(indexType vert1, indexType vert2)
    -> indexType {

  return this->edges.getEdgeNumberByVertexNumbers(vert1, vert2);
}

auto GeometryData::getFaceNumberByVerts(indexType vert1, indexType vert2,
                                        indexType vert3) -> indexType {

  return this->faces.getFaceNumberByVertexNumbers(vert1, vert2, vert3);
}

auto GeometryData::getxMax() -> Types::Vector3<prec> { return this->xMax; }

auto GeometryData::getxMin() -> Types::Vector3<prec> { return this->xMin; }

void GeometryData::sortReorientFacesPeriodicBC(
    PointerCollection &pointers, std::vector<indexType> &masterFaces,
    std::vector<indexType> &slaveFaces) {
  std::vector<Faces *> mFaces;
  std::vector<Faces *> sFaces;

  indexType numfaces = masterFaces.size();
  if (numfaces > 0) {
    if (numfaces != slaveFaces.size()) {
      throw std::runtime_error(
          "ERROR in GeometryData::sortReorientFacesPeriodicBC: Number of slave "
          "faces does not match number of master faces!");
    }

    mFaces.resize(numfaces);
    sFaces.resize(numfaces);
    for (indexType i = 0; i < numfaces; ++i) {
      mFaces[i] = this->getFace(masterFaces[i]);
      sFaces[i] = this->getFace(slaveFaces[i]);
    }

    // Make outward normal consistent
    Types::Vector3<prec> normal = mFaces[0]->getFaceNormal(pointers);
    for (auto i = mFaces.begin() + 1; i != mFaces.end(); ++i) {
      Types::Vector3<prec> tnormal = (*i)->getFaceNormal(pointers);
      if (abs(normal.dot(tnormal) - prec(1)) >
          std::numeric_limits<prec>::epsilon() * prec(1000))
        (*i)->flip();
    }
    for (auto i = sFaces.begin(); i != sFaces.end(); ++i) {
      Types::Vector3<prec> tnormal = (*i)->getFaceNormal(pointers);
      if (abs(normal.dot(tnormal) - prec(1)) >
          std::numeric_limits<prec>::epsilon() * prec(1000))
        (*i)->flip();
    }

    // compute local basis
    IntegrationPoint faceIp;
    faceIp.xi = prec(0);
    faceIp.eta = prec(0);
    Types::Matrix33<prec> localBasis;
    {
      Types::Vector3<prec> tangent;
      Types::Vector3<prec> tangent2;

      tangent = (mFaces[0]->getTangent_G1(pointers, faceIp)).normalized();
      tangent2 = normal.cross(tangent);
      localBasis.block<3, 1>(0, 0) = tangent;
      localBasis.block<3, 1>(0, 1) = tangent2;
      localBasis.block<3, 1>(0, 2) = normal;
    }

    // compute mean coordinate in local basis
    std::vector<Types::Vector3<prec>> masterMeanCoors;
    std::vector<Types::Vector3<prec>> slaveMeanCoors;
    masterMeanCoors.resize(numfaces);
    slaveMeanCoors.resize(numfaces);
    for (auto i = 0; i < numfaces; ++i) {
      masterMeanCoors[i] =
          localBasis.transpose() * mFaces[i]->computeMeanCoordinate(pointers);
      slaveMeanCoors[i] =
          localBasis.transpose() * sFaces[i]->computeMeanCoordinate(pointers);
    }

    // sort face with same mean coordinate x1 x2 in local basis.
    for (indexType i = 0; i < numfaces; ++i) {
      indexType j = i;
      bool search = true;
      while (search) {
        Types::Vector3<prec> diff = masterMeanCoors[i] - slaveMeanCoors[j];
        prec delt = abs(diff(0)) + abs(diff(1));
        if (delt < std::numeric_limits<prec>::epsilon() * prec(1000)) {
          std::swap(slaveMeanCoors[j], slaveMeanCoors[i]);
          std::swap(slaveFaces[j], slaveFaces[i]);
          search = false;
        }
        ++j;
        if (j == numfaces)
          search = false;
      }
    }

    // rotate face such that local coordinates match
    for (indexType i = 0; i < numfaces; ++i) {
      auto mFace = pointers.getGeometryData()->getFace(masterFaces[i]);
      auto sFace = pointers.getGeometryData()->getFace(slaveFaces[i]);

      auto mG1 = mFace->getTangent_G1(pointers, faceIp).normalized();
      auto sG1 = sFace->getTangent_G1(pointers, faceIp).normalized();

      prec dd = abs(mG1.dot(sG1) - prec(1));
      while (dd > std::numeric_limits<prec>::epsilon() * prec(1000)) {
        sFace->rotate(1);
        auto sG1 = sFace->getTangent_G1(pointers, faceIp).normalized();
        dd = abs(mG1.dot(sG1) - prec(1));
      }
    }

    // check edge orientation
    for (indexType i = 0; i < numfaces; ++i) {
      auto mFace = pointers.getGeometryData()->getFace(masterFaces[i]);
      auto sFace = pointers.getGeometryData()->getFace(slaveFaces[i]);

      indexType numEdges = mFace->getNumberOfEdges();
      for (indexType j = 0; j < numEdges; ++j) {
        auto mEdge = mFace->getEdge(pointers, j);
        auto sEdge = sFace->getEdge(pointers, j);

        auto mDir = mEdge->getA1Vector(pointers, faceIp);
        auto sDir = sEdge->getA1Vector(pointers, faceIp);

        prec dd = mDir.dot(sDir);
        if (dd < 0.99 || dd > 1.01) {
          sEdge->flip();
        }
      }
    }
  }
}

void GeometryData::checkUpdate() {

  for (auto face : this->faces) {
    face->checkUpdateElement(*this);
  }

  for (auto volume : this->volumes) {
    volume->checkUpdateElement(*this);
  }

  auto v1 = this->vertices.begin();
  this->xMin = v1->getCoordinates();
  this->xMax = this->xMin;

  for (auto v : this->vertices) {
    auto coor = v->getCoordinates();
    for (indexType i = 0; i < 3; ++i) {
      if (coor(i) > xMax(i))
        xMax(i) = coor(i);
      if (coor(i) < xMin(i))
        xMin(i) = coor(i);
    }
  }
}

} // namespace HierAMuS::Geometry
