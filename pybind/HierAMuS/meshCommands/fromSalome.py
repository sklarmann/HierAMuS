# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from HierAMuS.HierAMuSPyWrapper import HierAMuSPyFEM

class fromSalome:
    def __init__(self,meshcommands,_SMESH,_smesh):
        self.mesh = meshcommands
        self.SMESH = _SMESH
        self.smesh = _smesh

        self.geo = self.mesh.program.ptr.getGeometryData()
        pass

    def geomFromSalome(self,Mesh):
        nnIds = Mesh.GetNodesId()

        geoCMDs = self.mesh.getGeometryCommands()

        for i in nnIds:
            coor = Mesh.GetNodeXYZ(i)
            geoCMDs.addVertex(i,coor[0],coor[1],coor[2])


        elementNumbers = Mesh.GetElementsId()
        for i in elementNumbers:
            ElementType = Mesh.GetElementGeomType(i)

            if (ElementType == self.SMESH.Entity_Edge):  # Entity_Edge -> linearEdge
                verts = Mesh.GetElemNodes(i)
                geoCMDs.addLinearEdgeGeo(i,verts)

            elif (ElementType == self.SMESH.Entity_Quadrangle):  # Entity_Quadrangle -> quadrilateralFace
                verts = Mesh.GetElemNodes(i)
                #print(verts)
                edges=[]
                for j in range(3):
                    nn = [verts[j], verts[j+1]]
                    ee = Mesh.GetElementsByNodes(nn,self.SMESH.EDGE)
                    #print(nn)
                    #print(ee)
                    edges.append(ee[0])


                nn = [verts[3], verts[0]]
                ee = Mesh.GetElementsByNodes(nn,self.SMESH.EDGE)
                edges.append(ee[0])
                self.geo.requestNewGeometryObject(HierAMuSPyFEM.Geometry.GeometryTypes.LinearQuadrilateral,i)
                face = self.geo.getFace(i)
                face.setVerts(verts)
                face.setEdges(edges)

            elif (ElementType == self.SMESH.Entity_Hexa):
                nodeNums = Mesh.GetElemNodes(i)
                verts=[]
                for j in range(4):
                    verts.append(nodeNums[j+4])
                for j in range(4):
                    verts.append(nodeNums[j])

                edges =[]
        
                for j in range(3):
                    toAdd = Mesh.GetElementsByNodes(verts[j:j+2],self.SMESH.EDGE)
                    edges.append(toAdd[0])
                toAdd = Mesh.GetElementsByNodes([verts[3],verts[0]],self.SMESH.EDGE)
                edges.append(toAdd[0])

                for j in range(3):
                    toAdd = Mesh.GetElementsByNodes(verts[j+4:j+6],self.SMESH.EDGE)
                    edges.append(toAdd[0])
                toAdd = Mesh.GetElementsByNodes([verts[7],verts[4]],self.SMESH.EDGE)
                edges.append(toAdd[0])

                for j in range(4):
                    toAdd = Mesh.GetElementsByNodes([verts[0+j],verts[4+j]],self.SMESH.EDGE)
                    edges.append(toAdd[0])

                lll = [
                    [0,3,7,4],
                    [1,2,6,5],
                    [0,1,5,4],
                    [3,2,6,7],
                    [0,1,2,3],
                    [4,5,6,7]
                ]
                faces = []
                for ll in lll:
                    b=[]
                    for j in ll:
                        b.append(verts[j])
                    toAdd = Mesh.GetElementsByNodes(b,self.SMESH.FACE)
                    faces.append(toAdd[0])

                self.geo.requestNewGeometryObject(HierAMuSPyFEM.Geometry.GeometryTypes.LinearBrick,i)
                vol = self.geo.getVolume(i)
                vol.setVerts(verts)
                vol.setEdges(edges)
                vol.setFaces(faces)


    def addElementFromSalomeGroup(self,Group,MaterialNumber):
        elemCMDs = self.mesh.getElementCommands()
        Mesh = Group.GetMesh()
        elementList = Group.GetIDs()

        for i in elementList:
            elType = Mesh.GetElementGeomType(i)

            if(elType == self.SMESH.Entity_Edge):
                elemCMDs.addLinearEdge(MaterialNumber,i)

            elif(elType == self.SMESH.Entity_Quadrangle):
                elemCMDs.addFace(MaterialNumber,i)

            elif(elType == self.SMESH.Entity_Hexa):
                elemCMDs.addVolume(MaterialNumber,i)

    def SalomeTypeToFEMType(self,salomeType):
        if salomeType == self.SMESH.Entity_Edge:
            return HierAMuSPyFEM.Geometry.GeometryTypes.LinearEdge
        elif salomeType == self.SMESH.Entity_Quadrangle:
            return HierAMuSPyFEM.Geometry.GeometryTypes.LinearQuadrilateral
            
    def SalomeTypeToFEMGroupType(self,salomeType):
        if salomeType == self.SMESH.Entity_Edge:
            return HierAMuSPyFEM.Geometry.GeometryTypes.Edges
        elif salomeType == self.SMESH.Entity_Quadrangle:
            return HierAMuSPyFEM.Geometry.GeometryTypes.Faces

    def createInternalEdges(self,Mesh):
    

        # a = Mesh.GetElementsByType(SMESH.FACE)

        # for i in a:
        #     b = Mesh.GetElemNodes(i)
        #     b.append(b[0])
        #     for j in range(4):
        #         ee = [b[j], b[j+1]]
        #         Mesh.AddEdge(ee)

        a = Mesh.GetElementsByType(self.SMESH.VOLUME)
        for el in a:
            typ = Mesh.GetElementGeomType(el)
            if typ == self.SMESH.Entity_Hexa:
                nn = Mesh.GetElemNodes(el)

                ttemp = [[3,1,2,0], [4,5,6,7], [4,7,3,0], [5,6,2,1], [4,5,1,0], [7,6,2,3]]
                for temp in ttemp: 
                    b=[]
                    for i in temp:
                        b.append(nn[i])
                    Mesh.AddFace(b) 

        Mesh.MergeEqualElements()

        mgr = self.smesh.CreateFilterManager()
        functor = mgr.CreateMultiConnection2D()
        functor.SetMesh( Mesh.GetMesh() )
        for e in functor.GetValues():
            #if e.myNbConnects == 2:
            Mesh.AddEdge([e.myPnt1, e.myPnt2])

        Mesh.MergeEqualElements()