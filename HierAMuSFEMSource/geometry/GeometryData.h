// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once




#include "geometry/Vertex.h"
#include <forwarddeclaration.h>
#include <geometry/GeometryTypes.h>
#include <indexLists/rangedVector.h>
#include <stack>
#include <vector>
#include <unordered_map>

#include <geometry/datalists/vertexlist.h>
#include <geometry/datalists/edgelist.h>
#include <geometry/datalists/facelist.h>
#include <geometry/datalists/volumelist.h>

namespace HierAMuS::Geometry {

class GeometryData {
public:
  GeometryData();
  virtual ~GeometryData();
  auto requestNewVert() -> indexType;

  auto getVertex(indexType num) -> Vertex &;
  auto getEdge(indexType num) -> Edges &;
  auto getFace(indexType num) -> Faces *;
  auto getVolume(indexType num) -> Volumes *;
  auto getSpecial(indexType num) -> Special *;

  auto getVertexClosestTo(Types::Vector3<prec> point) -> Vertex &;

  auto requestNewGeometryObject(GeometryTypes type) -> indexType;
  auto requestNewGeometryObject(GeometryTypes type,
                                     indexType number) -> indexType;
  auto getGeometryElement(GeometryTypes type, indexType num) -> Base *;

  auto getNumberOfVertices() -> indexType {
    indexType ret = this->vertices.getNumberOfElements();
    return ret;
  };
  auto getNumberOfEdges() -> indexType {
    // indexType ret = this->edgeList.getNumberOfElements();
    indexType ret = this->edges.getNumberOfElements();
    return ret;
  };

  auto getLastVertexNumber() -> indexType { return this->vertices.lastElement(); };
  auto getNextVertexNumber() -> indexType { return this->vertices.lastElement() + 1; };

  void getGeometricElementInPlane(std::vector<prec> normal,
                                  std::vector<prec> point,
                                  GeometryTypes type,
                                  std::stack<Base *> &elems);

  
  auto getVerticesInPlane(PointerCollection& pointers, const Types::Vector3<prec>& normal, const Types::Vector3<prec>& point)
  -> std::vector<Geometry::Vertex *>;
  auto getEdgesInPlane(PointerCollection& pointers, const Types::Vector3<prec>& normal, const Types::Vector3<prec>& point)
  -> std::vector<Geometry::Edges *>;
  auto getFacesInPlane(PointerCollection& pointers, const Types::Vector3<prec>& normal, const Types::Vector3<prec>& point)
  -> std::vector<Geometry::Faces *>;

  void print(PointerCollection& pointers);

  void checkUpdate();

  auto getEdgeNumberByVerts(indexType vert1, indexType vert2) -> indexType;
  auto getFaceNumberByVerts(indexType vert1, indexType vert2, indexType vert3) -> indexType;

  auto getxMax() -> Types::Vector3<prec>;
  auto getxMin() -> Types::Vector3<prec>;


  /** @brief Sorts the slave face numbers such that they have the same local x1, x2 coordinates and are just shifted by x3.
    * Reorients the slave faces such that their local coordinate system matches the one of the master faces, this is required for 
    * periodic boundary conditions with hierarchical higher order shape functions.
    *
    *
    * @param[inout] masterFaces vector with the master face numbers.
    * @param[inout] slaveFaces vector with the slave face numbers.
    */
  void sortReorientFacesPeriodicBC(PointerCollection &pointers, std::vector<indexType> &masterFaces,
                                   std::vector<indexType> &slaveFaces);

private:
  //PointerCollection *pointers;

  vertexlist vertices;
  edgelist edges;
  facelist faces;
  volumelist volumes;

  Types::Vector3<prec> xMin, xMax;
  
  rangedVector<indexType, Special *> specialList;
};

} // namespace HierAMuS
