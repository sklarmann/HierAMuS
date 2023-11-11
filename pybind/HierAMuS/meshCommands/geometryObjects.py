# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from HierAMuS.HierAMuSPyWrapper import HierAMuSPyFEM
from HierAMuS.HierAMuSPyWrapper.HierAMuSPyFEM import Geometry

from timeit import default_timer as timer

class geometryObjects:
    def __init__(self,meshcmds):
        self.mesh = meshcmds
        self.geo = self.mesh.program.ptr.getGeometryData()
        self.ptr = self.mesh.program.ptr
        pass

    def vertexType(self):
        return Geometry.GeometryTypes.Vertex

    def edgeType(self):
        return Geometry.GeometryTypes.Edges

    def faceType(self):
        return Geometry.GeometryTypes.Faces

    def volumeType(self):
        return Geometry.GeometryTypes.Volumes

    def requestVertex(self,x,y,z):
        num = self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.Vertex);
        vert = self.geo.getVertex(num)
        vert.setCoordinates(x,y,z)
        return num

    def requestLinearEdgeGeo(self,vertexList):
        print("Deprecated function requestLinearEdgeGeo, wìll be renamed to requestLinearEdge")
        number = self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.LinearEdge)
        edge = self.geo.getEdge(number)
        edge.setVerts(self.geo,vertexList)
        return number

    def requestLinearEdge(self,vertexList):
        number = self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.LinearEdge)
        edge = self.geo.getEdge(number)
        edge.setVerts(self.geo,vertexList)
        return number
    
    
    def requestQuadraticEdge(self,vertexList):
        number = self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.QuadraticEdge)
        edge = self.geo.getEdge(number)
        edge.setVerts(self.geo,vertexList)
        return number

    def requestQuadrilateralFace(self,vertList,edgeList=None):
        print("Deprecated function requestQuadrilateralFace, will be renamed to requestLinearQuadrilateralFace")
        number = self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.LinearQuadrilateral)
        face = self.geo.getFace(number)
        face.setVerts(self.geo,vertList)
        if edgeList is not None:
            face.setEdges(edgeList)
        return number
    
    def requestLinearQuadrilateralFace(self,vertList,edgeList=None):
        number = self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.LinearQuadrilateral)
        face = self.geo.getFace(number)
        face.setVerts(self.geo,vertList)
        if edgeList is not None:
            face.setEdges(edgeList)
        return number
    
    def requestQuadraticQuadrilateralFace(self,vertList,edgeList=None):
        number = self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.QuadraticQuadrilateral)
        face = self.geo.getFace(number)
        face.setVerts(self.geo,vertList)
        if edgeList is not None:
            face.setEdges(edgeList)
        return number

    def requestTriangleFace(self,vertList,edgeList=None):
        print("Deprecated function requestTriangleFace, will be renamed to requestLinearTriangleFace")
        number = self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.LinearTriangle)
        face = self.geo.getFace(number)
        face.setVerts(self.geo,vertList)
        if edgeList is not None:
            face.setEdges(edgeList)
        return number
    
    def requestLinearTriangleFace(self,vertList,edgeList=None):
        number = self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.LinearTriangle)
        face = self.geo.getFace(number)
        face.setVerts(self.geo,vertList)
        if edgeList is not None:
            face.setEdges(edgeList)
        return number

    def requestLinearBrick(self,vertList,edgeList=None,faceList=None):
        number = self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.LinearBrick)
        vol = self.geo.getVolume(number)
        vol.setVerts(self.geo,vertList)
        if edgeList is not None:
            vol.setEdges(edgeList)
        if faceList is not None:
            vol.setFaces(faceList)
        return number

    def requestScaled2DFace(self,vertList,edgeList=None):
        number = self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.ScaledBoundary2D)
        face = self.geo.getFace(number)
        face.setVerts(self.geo,vertList)
        if edgeList is not None:
            face.setEdges(edgeList)
        return number

    def addVertex(self,num,x,y,z):
        self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.Vertex, num);
        vert = self.geo.getVertex(num)
        vert.setCoordinates(x,y,z)

    def addLinearEdgeGeo(self,number,vertexList):
        print("Deprecated function addLinearEdgeGeo, will be renamed to addLinearEdge")
        self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.LinearEdge,number)
        edge = self.geo.getEdge(number)
        edge.setVerts(self.geo,vertexList)
        
    def addLinearEdge(self,number,vertexList):
        self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.LinearEdge,number)
        edge = self.geo.getEdge(number)
        edge.setVerts(self.geo,vertexList)
        
    def addQuadraticEdge(self,number,vertexList):
        self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.QuadraticEdge,number)
        edge = self.geo.getEdge(number)
        edge.setVerts(self.geo,vertexList)
        
    def addQuadrilateralFace(self,number,vertList,edgeList = None):
        print("Deprecated function addQuadrilateralFace, will be renamed to addLinearQuadrilateralFace")
        self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.LinearQuadrilateral,number)
        face = self.geo.getFace(number)
        face.setVerts(self.geo,vertList)
        if edgeList is not None:
            face.setEdges(edgeList)
            
    def addLinearQuadrilateralFace(self,number,vertList,edgeList = None):
        self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.LinearQuadrilateral,number)
        face = self.geo.getFace(number)
        face.setVerts(self.geo,vertList)
        if edgeList is not None:
            face.setEdges(edgeList)
    
    def addQuadraticQuadrilateralFace(self,number,vertList,edgeList = None):
        self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.QuadraticQuadrilateral,number)
        face = self.geo.getFace(number)
        face.setVerts(self.geo,vertList)
        if edgeList is not None:
            face.setEdges(edgeList)

    def addTriangleFace(self,number,vertList,edgeList = None):
        print("Deprecated function addTriangleFace, will be renamed to addLinearTriangleFace")
        self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.LinearTriangle,number)
        face = self.geo.getFace(number)
        face.setVerts(self.geo,vertList)
        if edgeList is not None:
            face.setEdges(edgeList)
            
    def addLinearTriangleFace(self,number,vertList,edgeList = None):
        self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.LinearTriangle,number)
        face = self.geo.getFace(number)
        face.setVerts(self.geo,vertList)
        if edgeList is not None:
            face.setEdges(edgeList)

    def addLinearBrick(self,number,vertList,edgeList=None,faceList=None):
        self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.LinearBrick,number)
        vol = self.geo.getVolume(number)
        vol.setVerts(self.geo,vertList)
        if edgeList is not None:
            vol.setEdges(edgeList)
        if faceList is not None:
            vol.setFaces(faceList)

    def addScaled2DFace(self,number,vertList,edgeList=None):
        self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.ScaledBoundary2D,number)
        face = self.geo.getFace(number)
        face.setVerts(self.geo,vertList)
        if edgeList is not None:
            face.setEdges(edgeList)
        face.computeScalingCenter(self.ptr)

    def checkGeometry(self):
        start = timer()
        self.geo.checkUpdate(self.ptr.getEquationHandler())
        end = timer()
        self.mesh.program.outputLine(self.mesh.program.BasicLog(),self.mesh.program.BasicLog(),["Geometry check and update took:  ", end-start, "s\n"])

    def add2DBeamInterface(self,number,beamVertex,edgeList):
        self.geo.requestNewGeometryObject(self.ptr.getEquationHandler(),HierAMuSPyFEM.Geometry.GeometryTypes.BeamInterface2D,number)
        vol = self.geo.getSpecial(number)
        vol.setBeamVertex(beamVertex)
        vol.setEdges(edgeList)


