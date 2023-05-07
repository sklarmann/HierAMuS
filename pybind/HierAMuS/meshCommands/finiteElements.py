# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from HierAMuS.HierAMuSPyWrapper import HierAMuSPyFEM

class finiteElements:
    def __init__(self,meshcommands):
        self.mesh = meshcommands
        self.elements = self.mesh.program.ptr.getElementList()
        self.materials = self.mesh.program.ptr.getMaterialList()
        self.ptr = self.mesh.program.ptr
        #tempElem = HierAMuSPyFEM.FiniteElement.GenericFiniteElement()
        pass




    def addEdge(self,materialNumber,edge):
        if type(edge) == type([]):
            aedge = edge
        else:
            aedge = [edge]

        elem = self.elements.requestNewElement(self.ptr, HierAMuSPyFEM.FiniteElement.Elementtypes.Edge)
        elem.setEdges(aedge)
        elem.setMatrial(self.materials.getMaterial(materialNumber))

    def addFace(self,materialNumber,FaceNumber):
        if type(FaceNumber )==list:
            for i in FaceNumber:
                elem = self.elements.requestNewElement(self.ptr, HierAMuSPyFEM.FiniteElement.Elementtypes.Face)
                elem.setFace(i)
                elem.setMatrial(self.materials.getMaterial(materialNumber))


        else:
            elem = self.elements.requestNewElement(self.ptr, HierAMuSPyFEM.FiniteElement.Elementtypes.Face)
            elem.setFace(FaceNumber)
            elem.setMatrial(self.materials.getMaterial(materialNumber))
            
    
    def addFaceConstraint(self,materialNumber,FaceNumber,VertexNumber):
        if type(FaceNumber )==list:
            for i in FaceNumber:
                elem = self.elements.requestNewElement(self.ptr, HierAMuSPyFEM.FiniteElement.Elementtypes.FaceConstraint)
                elem.setFace(i)
                elem.setVerts([VertexNumber])
                elem.setMatrial(self.materials.getMaterial(materialNumber))
        else:
            elem = self.elements.requestNewElement(self.ptr, HierAMuSPyFEM.FiniteElement.Elementtypes.FaceConstraint)
            elem.setFace(FaceNumber)
            elem.setVerts([VertexNumber])
            elem.setMatrial(self.materials.getMaterial(materialNumber))


    def addVolume(self,materialNumber,VolumeNumber):
        if type(VolumeNumber )==list:
            for i in VolumeNumber:
                elem = self.elements.requestNewElement(self.ptr, HierAMuSPyFEM.FiniteElement.Elementtypes.Volume)
                elem.setVolume(i)
                elem.setMatrial(self.materials.getMaterial(materialNumber))


        else:
            elem = self.elements.requestNewElement(self.ptr, HierAMuSPyFEM.FiniteElement.Elementtypes.Volume)
            elem.setVolume(VolumeNumber)
            elem.setMatrial(self.materials.getMaterial(materialNumber))
            
    
    def addVolumeConstraint(self,materialNumber,VolumeNumber,VertexNumber):
        if type(VolumeNumber )==list:
            for i in VolumeNumber:
                elem = self.elements.requestNewElement(self.ptr, HierAMuSPyFEM.FiniteElement.Elementtypes.VolumeConstraint)
                elem.setVolume(i)
                elem.setVerts([VertexNumber])
                elem.setMatrial(self.materials.getMaterial(materialNumber))


        else:
            elem = self.elements.requestNewElement(self.ptr, HierAMuSPyFEM.FiniteElement.Elementtypes.VolumeConstraint)
            elem.setVolume(VolumeNumber)
            elem.setVerts([VertexNumber])
            elem.setMatrial(self.materials.getMaterial(materialNumber))



    def add2DInterfaceElement(self,materialNumber,specialNumber,vertNum,edgeList):
        geoData = self.mesh.program.ptr.getGeometryData()
        geoData.requestNewGeometryObject(HierAMuSPyFEM.Geometry.GeometryTypes.BeamInterface2D,specialNumber)
        elem = geoData.getSpecial(specialNumber)
        if type(vertNum) == list:
            elem.setBeamVertex(vertNum[0])
        else:
            elem.setBeamVertex(vertNum)

        elem.setEdges(edgeList)
        elem = self.elements.requestNewElement(self.ptr,HierAMuSPyFEM.FiniteElement.Elementtypes.beamInterfaceElement2D)
        if type(specialNumber) == list:
            elem.setSpecial(specialNumber)
        else:
            elem.setSpecial([specialNumber])
        elem.setMatrial(self.materials.getMaterial(materialNumber))


    def add3DInterfaceElement(self,materialNumber,specialNumber,vertNum,faceList,faceMaterialList):
        geoData = self.mesh.program.ptr.getGeometryData()
        geoData.requestNewGeometryObject(HierAMuSPyFEM.Geometry.GeometryTypes.BeamInterface3D,specialNumber)
        elem = geoData.getSpecial(specialNumber)
        if type(vertNum) == list:
            elem.setBeamVertex(vertNum[0])
        else:
            elem.setBeamVertex(vertNum)

        elem.setFaces(faceList)

        elem = self.elements.requestNewElement(self.ptr,HierAMuSPyFEM.FiniteElement.Elementtypes.beamInterfaceElement3D)
        if type(specialNumber) == list:
            elem.setSpecial(specialNumber)
        else:
            elem.setSpecial([specialNumber])
        elem.setMatrial(self.materials.getMaterial(materialNumber))
        elem.setFaces(faceList)
        elem.setMaterialPerSubElement(faceMaterialList)

        if type(vertNum) != list:
            vertNum = [vertNum]
        elem.setVerts(vertNum)

