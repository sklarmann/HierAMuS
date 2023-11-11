// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <cstddef>

#include <geometry/GeometryData.h>



#include "EquationHandler.h"

#include "geometry/GeometryBaseData.h"
#include "Edges/QuadraticEdge.h"
//#include "Faces/QuadraticQuadrilateral.h"
#include "Faces/ScaledBoundary2DData.h"
#include "geometry/Edges/EdgesRuntime.h"
#include "geometry/Faces/FacesRuntime.h"
#include "geometry/VertexRuntime.h"
#include "geometry/Volumes/VolumesRuntime.h"
#include <geometry/Edges/EdgesData.h>
#include <geometry/Edges/LinearEdgeData.h>
#include <geometry/Faces/FacesData.h>
#include <geometry/Faces/LinearQuadrilateralData.h>
#include <geometry/Faces/LinearTriangleData.h>
#include <geometry/VertexData.h>
#include <geometry/Volumes/LinearBrickData.h>
#include <geometry/Volumes/LinearPrismData.h>
#include <geometry/Volumes/VolumesData.h>

#include "geometry/Special/Special.h"
#include "geometry/Special/BeamInterface2D.h"

#include <math/Plane.h>

#include <algorithm>
#include <iterator>
#include <vector>

#include "geometry/datalists/GeometryList.h"
#include "geometry/datalists/GeometryListSingle.h"

namespace HierAMuS::Geometry {

GeometryData::GeometryData() {
  m_vertices = std::make_shared<GeometryListSingle<VertexData>>();
  m_edges = std::make_shared<GeometryList<EdgesData>>();
  m_faces = std::make_shared<GeometryList<FacesData>>();
  m_volumes = std::make_shared<GeometryList<VolumesData>>();
  m_special = std::make_shared<GeometryList<Special>>();
}

GeometryData::~GeometryData() {}

auto GeometryData::requestNewVert() -> indexType {
  indexType number = m_vertices->lastElement() + 1;
  m_vertices->add_element(number);
  return number;
}

auto GeometryData::getVertexData(indexType num) -> VertexData & {
  return m_vertices->get_element_reference(num);
}

auto GeometryData::getVertexRuntime(indexType num) -> VertexRuntime & {
  return *m_vertices_runtime.at(num);
}

auto GeometryData::getVertexRuntimePointer(indexType num) -> VertexRuntime * {
  return &(*m_vertices_runtime.at(num));
}

auto GeometryData::getEdgeData(indexType num) -> EdgesData & {
  auto &edge = m_edges->get_element_reference(num);
  edge.set_geometry_pointers(*this);
  return edge;
}

auto GeometryData::getEdgeRuntime(indexType num)
    -> std::shared_ptr<EdgesRuntime> {
  return m_edges_runtime.at(num);
}

auto GeometryData::getFaceData(indexType num) -> FacesData * {
  auto face = m_faces->get_element_pointer(num);
  face->set_geometry_pointers(*this);
  return face;
  ;
}

auto GeometryData::getFaceRuntime(indexType num)
    -> std::shared_ptr<FacesRuntime> {
  return m_faces_runtime.at(num);
}

auto GeometryData::getVolumeData(indexType num) -> VolumesData * {
  auto retVol = m_volumes->get_element_pointer(num);
  retVol->set_geometry_pointers(*this);
  return retVol;
}

auto GeometryData::getVolumeRuntime(indexType num)
    -> std::shared_ptr<VolumesRuntime> {
  return m_volumes_runtime.at(num);
}

auto GeometryData::getSpecial(indexType num) -> Special * {
  return m_special->get_element_pointer(num);
}

auto GeometryData::getVertexClosestTo(Types::Vector3<prec> point)
    -> VertexData & {
  prec minDist = std::numeric_limits<prec>::max();
  indexType minIndex = 0;

  for (auto &it : *m_vertices) {
    prec dist = (it.second.getCoordinates() - point).norm();
    if (dist < minDist) {
      minDist = dist;
      minIndex = it.first;
    }
  }
  return m_vertices->get_element_reference(minIndex);
}

auto GeometryData::requestNewGeometryObject(EquationHandler &eqHandler,
                                            GeometryTypes type) -> indexType {
  indexType ret = -1;
  switch (type) {
  case GeometryTypes::Vertex:
    ret = m_vertices->lastElement();
    break;
  case GeometryTypes::LinearEdge:
  case GeometryTypes::QuadraticEdge:
    ret = m_edges->lastElement();
    break;
  case GeometryTypes::LinearTriangle:
  case GeometryTypes::LinearQuadrilateral:
  case GeometryTypes::QuadraticQuadrilateral:
  case GeometryTypes::ScaledBoundary2D:
    ret = m_faces->lastElement();
    break;
  case GeometryTypes::LinearBrick:
  case GeometryTypes::LinearPrism:
    ret = m_volumes->lastElement();
    break;
  case GeometryTypes::BeamInterface2D:
  //case GeometryTypes::BeamInterface3D:
    ret = m_special->lastElement();
    break;
  default:
    throw std::invalid_argument("Selected Geometry type to add does not exist");
    break;
  }
  ++ret;
  this->requestNewGeometryObject(eqHandler, type, ret);
  return ret;
}

auto GeometryData::requestNewGeometryObject(EquationHandler &eqHandler,
                                            GeometryTypes type,
                                            indexType number) -> indexType {
  switch (type) {
  case GeometryTypes::Vertex:
    m_vertices->add_element(number);
    m_vertices->get_element_reference(number).setNodeSetManager(
        eqHandler.getNewNodeSetManager());
    break;
  case GeometryTypes::LinearEdge:
    m_edges->add_element<LinearEdgeData>(number);
    m_edges->get_element_reference(number).setNodeSetManager(
        eqHandler.getNewNodeSetManager());
    break;
  case GeometryTypes::QuadraticEdge:
    m_edges->add_element<QuadraticEdge>(number);
    m_edges->get_element_reference(number).setNodeSetManager(
        eqHandler.getNewNodeSetManager());
    break;
  case GeometryTypes::LinearTriangle:
    m_faces->add_element<LinearTriangleData>(number);
    m_faces->get_element_reference(number).setNodeSetManager(
        eqHandler.getNewNodeSetManager());
    break;
  case GeometryTypes::LinearQuadrilateral:
    m_faces->add_element<LinearQuadrilateralData>(number);
    m_faces->get_element_reference(number).setNodeSetManager(
        eqHandler.getNewNodeSetManager());
    break;
  //case GeometryTypes::QuadraticQuadrilateral:
  //  m_faces->add_element<QuadraticQuadrilateral>(number);
  //  m_faces->get_element_reference(number).setNodeSetManager(
  //      pointers.getEquationHandler()->getNewNodeSetManager());
  //  break;
  case GeometryTypes::ScaledBoundary2D:
    m_faces->add_element<ScaledBoundary2DData>(number);
    m_faces->get_element_reference(number).setNodeSetManager(
        eqHandler.getNewNodeSetManager());
    break;
  case GeometryTypes::LinearBrick:
    m_volumes->add_element<LinearBrickData>(number);
    m_volumes->get_element_reference(number).setNodeSetManager(
        eqHandler.getNewNodeSetManager());
    break;
  case GeometryTypes::LinearPrism:
    m_volumes->add_element<LinearPrismData>(number);
    m_volumes->get_element_reference(number).setNodeSetManager(
        eqHandler.getNewNodeSetManager());
    break;
  case GeometryTypes::BeamInterface2D:
    m_special->add_element<BeamInterface2D>(number);
    m_special->get_element_reference(number).setNodeSetManager(
        eqHandler.getNewNodeSetManager());
    break;
  /** case GeometryTypes::BeamInterface3D:
    m_special->add_element<BeamInterface3D>(number);
    m_special->get_element_reference(number).setNodeSetManager(
        eqHandler.getNewNodeSetManager());
    break;**/
  default:
    // TODO throw exception
    throw std::runtime_error("In GeometryData->requstNewGeometryObject(), "
                             "tried to add unknown geometry type!");
    break;
  }
  return number;
}

auto GeometryData::getGeometryElement(GeometryTypes type, indexType num)
    -> GeometryBaseData * {

  switch (type) {
  case GeometryTypes::Vertex:
    return m_vertices->get_element_pointer(num);
    break;
  case GeometryTypes::Edges:
    return m_edges->get_element_pointer(num);
    break;
  case GeometryTypes::Faces:
    return m_faces->get_element_pointer(num);
    break;
  case GeometryTypes::Volumes:
    return m_volumes->get_element_pointer(num);
    break;
  case GeometryTypes::Special:
    return m_special->get_element_pointer(num);
    break;
  default:
    // TODO throw exception
    std::cout << "getGeometryElement Error" << std::endl;
  }
  return nullptr;
}

void GeometryData::getGeometricElementInPlane(
    std::vector<prec> normal, std::vector<prec> point, GeometryTypes type,
    std::stack<GeometryBaseData *> &elems) {

  Plane plane;
  plane.set(normal, point);
  std::vector<indexType> vertNums;

  for (auto vertice : *m_vertices) {
    if (plane.inPlane(vertice.second.getCoordinates())) {
      vertNums.push_back(vertice.second.getId());
    }
  }

  switch (type) {
  case GeometryTypes::Vertex: {
    for (auto &vertNum : vertNums) {
      // elems.push(this->vertList[*it]);
      elems.push(m_vertices->get_element_pointer(vertNum));
    }
  } break;
  case GeometryTypes::Edges: {
    if (!vertNums.empty()) {
      std::vector<indexType> edgeVerts;
      for (auto &edge : *m_edges) {
        bool add = true;
        edgeVerts.clear();
        edgeVerts = edge.second->getVertexNumbers();
        for (auto &edgeVert : edgeVerts) {
          if (std::find(vertNums.begin(), vertNums.end(), edgeVert) ==
              vertNums.end()) {
            add = false;
          }
        }
        if (add) {
          elems.push(&(*edge.second));
        }
      }
    }
    break;
  case GeometryTypes::Faces: {
    if (!vertNums.empty()) {
      std::vector<indexType> faceVerts;
      for (auto &face : *m_faces) {
        faceVerts.clear();
        faceVerts = face.second->getVertexNumbers();
        bool add = true;
        for (auto &faceVert : faceVerts) {
          if (std::find(vertNums.begin(), vertNums.end(), faceVert) ==
              vertNums.end()) {
            add = false;
          }
        }
        if (add) {
          elems.push(&(*face.second));
        }
      }
    }
  }
  } break;
  default:
    break;
  }
}

auto GeometryData::getVerticesInPlane(const Types::Vector3<prec> &normal,
                                      const Types::Vector3<prec> &point)
    -> std::vector<Geometry::VertexData *> {
  Plane plane;
  plane.set(normal, point);

  std::vector<VertexData *> verts;
  for (auto &v : *m_vertices) {
    if (plane.inPlane(v.second.getCoordinates()))
      verts.push_back(&(v.second));
  }
  return verts;
}

auto GeometryData::getEdgesInPlane(const Types::Vector3<prec> &normal,
                                   const Types::Vector3<prec> &point)
    -> std::vector<Geometry::EdgesData *> {
  Plane plane;
  plane.set(normal, point);

  auto verts = this->getVerticesInPlane(normal, point);

  std::map<indexType, EdgesData *> idEdgeMap;
  IntegrationPoint ip;
  ip.xi = 0;
  for (auto v : verts) {
    auto pedges = v->getConnectedEdges();
    for (auto edgenum : pedges) {
      auto &edge = this->getEdgeData(edgenum);
      if (plane.inPlane(edge.getCoordinates(ip)))
        idEdgeMap[edge.getId()] = &edge;
    }
  }

  std::vector<Geometry::EdgesData *> edgeList;
  for (auto j : idEdgeMap) {
    edgeList.push_back(j.second);
    j.second->set_geometry_pointers(*this);
  }

  return edgeList;
}

auto GeometryData::getFacesInPlane(const Types::Vector3<prec> &normal,
                                   const Types::Vector3<prec> &point)
    -> std::vector<Geometry::FacesData *> {
  Plane plane;
  plane.set(normal, point);

  auto verts = this->getVerticesInPlane(normal, point);

  std::map<indexType, FacesData *> idFaceMap;
  IntegrationPoint ip;
  ip.xi = 0;
  ip.eta = 0;
  for (auto v : verts) {
    auto pfaces = v->getConnectedFaces();
    for (auto facenum : pfaces) {
      auto face = this->getFaceData(facenum);
      if (plane.inPlane(face->getCoordinates(ip)))
        idFaceMap[face->getId()] = face;
    }
  }

  std::vector<Geometry::FacesData *> faceList;
  for (auto j : idFaceMap) {
    faceList.push_back(j.second);
  }

  return faceList;
}

void GeometryData::print(spdlog::logger &Logger) {


  Logger.info("\nGeometry Data informations:\n"
              "   Number of vertices:   {:>12}\n"
              "   Number of edges:      {:>12}\n"
              "   Number of faces:      {:>12}\n"
              "   Number of volumes:    {:>12}",
              m_vertices->getNumberOfElements(), m_edges->getNumberOfElements(),
              m_faces->getNumberOfElements(), m_volumes->getNumberOfElements());

  for (auto &vertice : *m_vertices) {
    vertice.second.print(Logger);
    Logger.debug("{:-<100}", "");
  }
  for (auto &edge : *m_edges) {
    edge.second->print(Logger);
    Logger.debug("{:-<100}", "");
  }
  for (auto &face : *m_faces) {
    face.second->print(Logger);
    Logger.debug("{:-<100}", "");
  }
  for (auto &vol : *m_volumes) {
    vol.second->print(Logger);
    Logger.debug("{:-<100}", "");
  }
}

auto GeometryData::getEdgeNumberByVerts(indexType vert1, indexType vert2)
    -> indexType {
  for (auto &edge : *m_edges) {
    if (edge.second->hasVertices(vert1, vert2))
      return edge.second->getId();
  }
  return -1;
}

auto GeometryData::getFaceNumberByVerts(indexType vert1, indexType vert2,
                                        indexType vert3) -> indexType {

  for (auto &face : *m_faces) {
    if (face.second->hasVertices(vert1, vert2, vert3)) {
      return face.second->getId();
    }
  }
  return -1;
}

auto GeometryData::getxMax() -> Types::Vector3<prec> { return this->xMax; }

auto GeometryData::getxMin() -> Types::Vector3<prec> { return this->xMin; }

void GeometryData::sortReorientFacesPeriodicBC(
    std::vector<indexType> &masterFaces,
    std::vector<indexType> &slaveFaces) {
  std::vector<FacesData *> mFaces;
  std::vector<FacesData *> sFaces;

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
      mFaces[i] = this->getFaceData(masterFaces[i]);
      sFaces[i] = this->getFaceData(slaveFaces[i]);
    }

    // Make outward normal consistent
    Types::Vector3<prec> normal = mFaces[0]->getFaceNormal();
    for (auto i = mFaces.begin() + 1; i != mFaces.end(); ++i) {
      Types::Vector3<prec> tnormal = (*i)->getFaceNormal();
      if (abs(normal.dot(tnormal) - prec(1)) >
          std::numeric_limits<prec>::epsilon() * prec(1000)) {
        (*i)->flip();
        auto runtime = this->getFaceRuntime((*i)->getId());
        runtime->flip();
      }
    }
    for (auto i = sFaces.begin(); i != sFaces.end(); ++i) {
      Types::Vector3<prec> tnormal = (*i)->getFaceNormal();
      if (abs(normal.dot(tnormal) - prec(1)) >
          std::numeric_limits<prec>::epsilon() * prec(1000)) {
        (*i)->flip();
        auto runtime = this->getFaceRuntime((*i)->getId());
        runtime->flip();
      }
    }

    // compute local basis
    IntegrationPoint faceIp;
    faceIp.xi = prec(0);
    faceIp.eta = prec(0);
    Types::Matrix33<prec> localBasis;
    {
      Types::Vector3<prec> tangent;
      Types::Vector3<prec> tangent2;

      tangent = (mFaces[0]->getTangent_G1(faceIp)).normalized();
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
          localBasis.transpose() * mFaces[i]->computeMeanCoordinate();
      slaveMeanCoors[i] =
          localBasis.transpose() * sFaces[i]->computeMeanCoordinate();
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
      auto mFace = this->getFaceData(masterFaces[i]);
      auto sFace = this->getFaceData(slaveFaces[i]);

      auto mG1 = mFace->getTangent_G1(faceIp).normalized();
      auto sG1 = sFace->getTangent_G1(faceIp).normalized();

      prec dd = abs(mG1.dot(sG1) - prec(1));
      while (dd > std::numeric_limits<prec>::epsilon() * prec(1000)) {
        sFace->rotate(1);
        this->getFaceRuntime(sFace->getId())->rotate(1);
        auto sG1 = sFace->getTangent_G1(faceIp).normalized();
        dd = abs(mG1.dot(sG1) - prec(1));
      }
    }

    // check edge orientation
    for (indexType i = 0; i < numfaces; ++i) {
      auto mFace = this->getFaceData(masterFaces[i]);
      auto sFace = this->getFaceData(slaveFaces[i]);

      indexType numEdges = mFace->getNumberOfEdges();
      for (indexType j = 0; j < numEdges; ++j) {
        auto mEdge = mFace->getEdge(j);
        auto sEdge = sFace->getEdge(j);

        auto mDir = mEdge->getA1Vector(faceIp);
        auto sDir = sEdge->getA1Vector(faceIp);

        prec dd = mDir.dot(sDir);
        if (dd < 0.99 || dd > 1.01) {
          sEdge->flip();
          auto runtime = this->getEdgeRuntime(sEdge->getId());
          runtime->flip();
        }
      }
    }
  }
}

void GeometryData::createRuntimeObjects() {
  for (auto &it : *m_vertices) {
    m_vertices_runtime.try_emplace(it.first, std::make_shared<VertexRuntime>(it.second));
  }
  for (auto &it : *m_edges) {
    m_edges_runtime.try_emplace(it.first,
                                it.second->getRuntimeObject(*this));
  }
  for (auto &it : *m_faces) {
    m_faces_runtime.try_emplace(it.first,
                                it.second->getRuntimeObject(*this));
  }
  for (auto &it : *m_volumes) {
    m_volumes_runtime.try_emplace(it.first,
                                  it.second->getRuntimeObject(*this));
  }
}

void GeometryData::updateRuntimeObjectEquations() {
  for (auto &it : *m_vertices) {
    m_vertices_runtime.at(it.first)->updateNodes(it.second.getNodeListMap());
  }
  for (auto &it : *m_edges) {
    m_edges_runtime.at(it.first)->updateNodes(it.second->getNodeListMap());
  }
  for (auto &it : *m_faces) {
    m_faces_runtime.at(it.first)->updateNodes(it.second->getNodeListMap());
  }
  for (auto &it : *m_volumes) {
    m_volumes_runtime.at(it.first)->updateNodes(it.second->getNodeListMap());
  }
}

void GeometryData::checkUpdate(EquationHandler &eqHandler) {

  for (auto &face : *m_faces) {
    face.second->checkUpdateElement(eqHandler, *this);
  }

  for (auto &volume : *m_volumes) {
    volume.second->checkUpdateElement(eqHandler, *this);
  }

  auto v1 = m_vertices->begin();
  this->xMin = v1->second.getCoordinates();
  this->xMax = this->xMin;

  for (auto &v : *m_vertices) {
    auto coor = v.second.getCoordinates();
    for (indexType i = 0; i < 3; ++i) {
      if (coor(i) > xMax(i))
        xMax(i) = coor(i);
      if (coor(i) < xMin(i))
        xMin(i) = coor(i);
    }
  }
}

auto GeometryData::getNextVertexNumber() -> indexType {
  return m_vertices->lastElement() + 1;
}

auto GeometryData::getLastVertexNumber() -> indexType {
  return m_vertices->lastElement();
}

auto GeometryData::getNumberOfEdges() -> indexType {
  // indexType ret = this->edgeList.getNumberOfElements();
  indexType ret = m_edges->getNumberOfElements();
  return ret;
}

auto GeometryData::getNumberOfVertices() -> indexType {
  indexType ret = m_vertices->getNumberOfElements();
  return ret;
}

} // namespace HierAMuS::Geometry
