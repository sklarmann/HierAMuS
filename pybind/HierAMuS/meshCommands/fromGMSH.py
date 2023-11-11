# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from HierAMuS.HierAMuSPyWrapper import HierAMuSPyFEM
import numpy as np
from timeit import default_timer as timer

class fromGMSh:
    def __init__(self,meshcmd) -> None:
        self.mesh = meshcmd
        pass


    def addGeomFromGmsh(self,gmsh):
        start = timer()
        self.addVertices(gmsh)

        eltypes, eltags, elnodes = gmsh.model.mesh.getElements()

        for i in range(len(eltypes)):
            eltype = eltypes[i]
            eltag = eltags[i]
            elnode = elnodes[i]

            if eltype == gmsh.model.mesh.getElementType("Line", 1):
                for j in range(len(eltag)):
                    self.mesh.getGeometryCommands().addLinearEdge(eltag[j], elnode[j*2:j*2+2])
            elif eltype == gmsh.model.mesh.getElementType("Line", 2):
                for j in range(len(eltag)):
                    telnode = [ elnode[j*3], elnode[j*3+1], elnode[j*3+2]]
                    self.mesh.getGeometryCommands().addQuadraticEdge(eltag[j],telnode)
            elif eltype == gmsh.model.mesh.getElementType("Quadrangle", 1):
                for j in range(len(eltag)):
                    self.mesh.getGeometryCommands().addLinearQuadrilateralFace(eltag[j], elnode[j*4:j*4+4])
            elif eltype == gmsh.model.mesh.getElementType("Quadrangle", 2):
                for j in range(len(eltag)):
                    self.mesh.getGeometryCommands().addQuadraticQuadrilateralFace(eltag[j], elnode[j*9:j*9+9])
            elif eltype == gmsh.model.mesh.getElementType("Triangle", 1):
                for j in range(len(eltag)):
                    self.mesh.getGeometryCommands().addLinearTriangleFace(eltag[j], elnode[j*3:j*3+3])
            elif eltype == gmsh.model.mesh.getElementType("Hexahedron", 1):
                for j in range(len(eltag)):
                    self.mesh.getGeometryCommands().addLinearBrick(eltag[j], elnode[j*8:j*8+8])


        end = timer()
        self.mesh.program.outputLine(self.mesh.program.BasicLog(),self.mesh.program.BasicLog(),["Time to add geometry from Gmsh: ", end-start, "s\n"])
        #print("Time to add geometry from Gmsh: ", end-start, "s\n")




    def addVertices(self,gmsh):
        nodetags, coords, parCoor = gmsh.model.mesh.getNodes()
        for i in range(len(nodetags)):
            self.mesh.getGeometryCommands().addVertex(num=nodetags[i], x=coords[i*3+0], y=coords[i*3+1], z=coords[i*3+2])


    def addLinearEdgesGeo(self,gmsh,order,lineTags=None):

        ltags = []
        if type(lineTags) == list:
            ltags = lineTags
        elif lineTags == None:
            ltags.append(-1)
        else:
            ltags.append(lineTags)

        geo = self.mesh.getGeometryCommands()
        for tag in ltags:
            eltype = gmsh.model.mesh.getElementType("Line", order)
            edgeTags, edgeNodes = gmsh.model.mesh.getElementsByType(eltype,tag=tag)
            numEdgeNodes = order+1

            if order==1:
                for i in range(len(edgeTags)):
                    enum = edgeTags[i]
                    en = edgeNodes[i*numEdgeNodes:i*numEdgeNodes+numEdgeNodes]
                    geo.addLinearEdgeGeo(enum, en)
                    
    def addLineElements(self,gmsh,lineTags,material):
        ltags = []
        if type(lineTags) == list:
            ltags = lineTags
        else:
            ltags.append(lineTags)
            
        
        elems = self.mesh.getElementCommands()
        for tag in ltags:
            [et, enum, nt] = gmsh.model.mesh.getElements(1,tag)
            for j in enum[0]:
                elems.addEdge(materialNumber=material,edge=j)
        

    def addFaceElements(self,gmsh,faceTags,material):
        ltags = []
        if type(faceTags) == list:
            ltags = faceTags
        else:
            ltags.append(faceTags)

        elems = self.mesh.getElementCommands()
        for tag in ltags:
            [et, enum, nt] = gmsh.model.mesh.getElements(2,tag)
            for j in enum[0]:
                elems.addFace(materialNumber=material, FaceNumber=j)
    
    
    def addFaceElementsPhysGroup(self,gmsh,groupName:str,material:int):
        tags = self.getTagsByPhysicalGroupName(gmsh,groupName,2)
        if tags:
            self.addFaceElements(gmsh,tags,material)          
    

    def addVolumeElements(self,gmsh,volTags,material:int):
        ltags = []
        if type(volTags) == list:
            ltags = volTags
        else:
            ltags.append(volTags)
        
        elems = self.mesh.getElementCommands()
        for tag in ltags:
            [et, enum, nt] = gmsh.model.mesh.getElements(3,tag)
            for j in enum[0]:
                elems.addVolume(materialNumber=material, VolumeNumber=j)
    
    def addVolumeElementsPhysGroup(self,gmsh,groupName:str,material:int):
        tags = self.getTagsByPhysicalGroupName(gmsh,groupName,3)
        if tags:
            self.addVolumeElements(gmsh,tags,material)          
        
    def addFaceConstraint(self,gmsh,faceTags,vertexTag,material):
        ltags = []
        if type(faceTags) == list:
            ltags = faceTags
        else:
            ltags.append(faceTags)
            
        [coord, pcoord, dim, nt] = gmsh.model.mesh.getNode(vertexTag)
        nn = nt
        elems = self.mesh.getElementCommands()
        for tag in ltags:
            [et, enum, nt] = gmsh.model.mesh.getElements(2,tag)
            for j in enum[0]:
                elems.addFaceConstraint(materialNumber=material,FaceNumber=j,VertexNumber=nn)
        
    def addVolumeConstraint(self,gmsh,volumeTags,vertexTag,material):
        ltags = []
        if type(volumeTags) == list:
            ltags = volumeTags
        else:
            ltags.append(volumeTags)
            
        [coord, pcoord, dim, nt] = gmsh.model.mesh.getNode(vertexTag)
        nn = nt
        elems = self.mesh.getElementCommands()
        for tag in ltags:
            [et, enum, nt] = gmsh.model.mesh.getElements(3,tag)
            for j in enum[0]:
                elems.addVolumeConstraint(materialNumber=material,VolumeNumber=j,VertexNumber=nn)
        


    def addQuadrilateralFaceGeo(self,gmsh,order,faceTags=None):
        ltags = []
        if type(faceTags) == list:
            ltags = faceTags
        elif faceTags == None:
            ltags.append(-1)
        else:
            ltags.append(faceTags)

        geo = self.mesh.getGeometryCommands()

        for tag in ltags:
            eltype = gmsh.model.mesh.getElementType("Quadrangle", order)
            facetags, facenodes = gmsh.model.mesh.getElementsByType(eltype,tag=tag)
            numFaceNodes = (order+1)**2

            if order == 1:
                for i in range(len(facetags)):
                    fnum = facetags[i]
                    fn = facenodes[i*numFaceNodes:(i+1)*numFaceNodes]
                    geo.addQuadrilateralFace(fnum, fn)


    def addQuadrilateralFiniteElements(self,gmsh,order,faceTags,material):
        ltags = []
        if type(faceTags) == list:
            ltags = faceTags
        else:
            ltags.append(faceTags)


        elems = self.mesh.getElementCommands()
        eltype = gmsh.model.mesh.getElementType("Quadrangle", order)

        for tag in ltags:
            facetags, facenodes = gmsh.model.mesh.getElementsByType(eltype,tag=tag)

            for i in range(len(facetags)):
                fnum = facetags[i]
                elems.addFace(materialNumber=material, FaceNumber=fnum)

    def addTriangleFiniteElements(self,gmsh,order,faceTags,material):
        ltags = []
        if type(faceTags) == list:
            ltags = faceTags
        else:
            ltags.append(faceTags)


        elems = self.mesh.getElementCommands()
        eltype = gmsh.model.mesh.getElementType("Triangle", order)

        for tag in ltags:
            facetags, facenodes = gmsh.model.mesh.getElementsByType(eltype,tag=tag)

            for i in range(len(facetags)):
                fnum = facetags[i]
                elems.addFace(materialNumber=material, FaceNumber=fnum)

    def addBrickVolumeGeo(self,gmsh,order,volumeTags=None):
        ltags = []
        if type(volumeTags) == list:
            ltags = volumeTags
        elif volumeTags == None:
            ltags.append(-1)
        else:
            ltags.append(volumeTags)

        geo = self.mesh.getGeometryCommands()

        eltype = gmsh.model.mesh.getElementType("Hexahedron", order)
        nodesPerElement = (order+1)**3

        for tag in ltags:
            voltags, volnodes = gmsh.model.mesh.getElementsByType(eltype,tag=tag)
            for i in range(len(voltags)):
                voln = volnodes[nodesPerElement*i:nodesPerElement*i+nodesPerElement]
                vt = voltags[i]
                geo.addBrickFiniteElementGeo(vt, voln)

    def addBrickVolumeFiniteElements(self,gmsh,order,volumeTags,material):
        ltags = []
        if type(volumeTags) == list:
            ltags = volumeTags
        else:
            ltags.append(volumeTags)

        geo = self.mesh.getGeometryCommands()

        eltype = gmsh.model.mesh.getElementType("Hexahedron", order)
        elems = self.mesh.getElementCommands()

        for tag in ltags:
            voltags, volnodes = gmsh.model.mesh.getElementsByType(eltype,tag=tag)
            for i in range(len(voltags)):
                vt = voltags[i]
                elems.addVolume(materialNumber=material, VolumeNumber=vt)







    def addQuadrilateralGeoElements(self,gmsh,order=1):
        facetags, facenodes = gmsh.model.mesh.getAllFaces(4)
        #eltyp = gmsh.model.mesh.getElementType("Quadrangle", 1)
        #facetags, facenodes = gmsh.model.mesh.getAllFaces(eltyp)

        for i in range(len(facetags)):
            fn = facenodes[i*4:i*4+4]
            el = fn[0:2]
            el = np.append(el,fn[1:3])
            el = np.append(el,fn[2:4])
            el = np.append(el,[fn[3],fn[0]])
            et,eo = gmsh.model.mesh.getEdges(el)
            self.mesh.getGeometryCommands().addQuadrilateralFace(number=facetags[i], vertList=fn, edgeList=et)




    def addFiniteElementsFromTag(self,gmsh,eltype,tag,material):
        et,nt = gmsh.model.mesh.getElementsByType(elementType=eltype,tag=tag)

        if eltype == gmsh.model.mesh().getElementType("Line", 1):
            et = gmsh.model.mesh.getEdges(nt)
            self.mesh.getElementCommands().addLinearEdge(material, et.tolist())
        elif gmsh.model.mesh().getElementType("Quadrangle", 1):
            et,fo = gmsh.model.mesh.getFaces(4,nt)
            self.mesh.getElementCommands().addFace(material, et.tolist())
        else:
            pass



    def getVertexNumbers(self,gmsh,dim,tag):
        nt,coord,paramCoord=gmsh.model.mesh.getNodes(dim,tag,includeBoundary=True)
        return nt


    def getEdgeNumbers(self,gmsh,tag,order=1):
        eltyp = gmsh.model.mesh.getElementType("Line", order)
        et, nt = gmsh.model.mesh.getElementsByType(eltyp,tag=tag)
        return et.tolist()

    def getFaceNumbers(self,gmsh,tagin,ftype,order):
        tag = []
        if type(tagin) != list:
            tag.append(tagin)
        else:
            tag = tagin

        eltyp = 0
        if ftype == 3:
            eltyp = gmsh.model.mesh.getElementType("Triangle", order)
        elif ftype == 4:
            eltyp = gmsh.model.mesh.getElementType("Quadrangle", order)
        else:
            raise Exception("Unknown face type")

        facelist = []
        for i in tag:
            et, nt = gmsh.model.mesh.getElementsByType(eltyp,tag=i)
            facelist.extend(et.tolist())
        return facelist
    
    def getTagsByPhysicalGroupName(self,gmsh,groupName:str,dim:int):
        groups = gmsh.model.getPhysicalGroups(dim)
        
        for i in groups:
            cname = gmsh.model.getPhysicalName(dim,i[1])
            if cname == groupName:
                return gmsh.model.getEntitiesForPhysicalGroup(dim,i[1]).tolist()

