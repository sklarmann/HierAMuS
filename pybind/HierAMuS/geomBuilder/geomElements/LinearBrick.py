# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from HierAMuS.geomBuilder.geomElements.Vertex import Vertex
from HierAMuS.geomBuilder.geomElements.LinearLine import Line

class LinearBrick:
    def __init__(self,fesys):
        self.fem = fesys
        self.geoId = -1
        self.vertList = []
        self.hEdgeList = []
        self.vEdgeList = []
        self.quadList = []

    def plotQuadLine(self,quadnum,ax):
        q = self.fem.ptr.getGeometryData().getFace(quadnum)
        vnums = q.getVertNums()

        verts = []
        for i in vnums:
            verts.append(self.fem.ptr.getGeometryData().getVertex(i))

        [x1,y1,z1] = verts[0].getCoordinates()
        [x2,y2,z2] = verts[1].getCoordinates()
        [x3,y3,z3] = verts[2].getCoordinates()
        [x4,y4,z4] = verts[3].getCoordinates()

        x = [x1,x2,x3,x4]
        y = [y1,y2,y3,y4]
        z = [z1,z2,z3,z4]
        ax.plot(x,y,z,color='b')
        x = [x1,x4]
        y = [y1,y4]
        z = [z1,z4]
        ax.plot(x,y,z,color='b')
        

    def plotQuadLineLine(self,quadnum,ax):
        q = self.fem.ptr.getGeometryData().getFace(quadnum)
        edgeNums = q.getEdgeNums()

        for i in edgeNums:
            edge = self.fem.ptr.getGeometryData().getEdge(i)
            vnums = edge.getVertNums()
            verts = []
            for j in vnums:
                verts.append(self.fem.ptr.getGeometryData().getVertex(j))

            c1 = verts[0].getCoordinates()
            c2 = verts[1].getCoordinates()
            [x,y,z] = [[c1[0],c2[0]],[c1[1],c2[1]],[c1[2],c2[2]]]
            ax.plot(x,y,z,color='b')


    def setElementSize(self,esx,esy,esz):
        Lx1 = self.edges[0].getLength()
        Lx2 = self.edges[2].getLength()
        Lx3 = self.edges[8].getLength()
        Lx4 = self.edges[10].getLength()

        Ly1 = self.edges[1].getLength()
        Ly2 = self.edges[3].getLength()
        Ly3 = self.edges[9].getLength()
        Ly4 = self.edges[11].getLength()

        Lz1 = self.edges[4].getLength()
        Lz2 = self.edges[5].getLength()
        Lz3 = self.edges[6].getLength()
        Lz4 = self.edges[7].getLength()
        
        Lax = (Lx1+Lx2+Lx3+Lx4)/4
        Lay = (Ly1+Ly2+Ly3+Ly4)/4
        Laz = (Lz1+Lz2+Lz3+Lz4)/4
        
        ndivx = round(Lax/esx)
        ndivy = round(Lay/esy)
        ndivz = round(Laz/esz)
        
        self.setDivision(ndivx,ndivy,ndivz)


    def setDivision(self,nx,ny,nz):
        self.nx = nx
        self.ny = ny
        self.nz = nz

        self.edges[0].setNumElements(nx)
        self.edges[1].setNumElements(ny)
        self.edges[2].setNumElements(nx)
        self.edges[3].setNumElements(ny)

        self.edges[4].setNumElements(nz)
        self.edges[5].setNumElements(nz)
        self.edges[6].setNumElements(nz)
        self.edges[7].setNumElements(nz)

        self.edges[8].setNumElements(nx)
        self.edges[9].setNumElements(ny)
        self.edges[10].setNumElements(nx)
        self.edges[11].setNumElements(ny)


    def setVerts(self,vertList:list):
        self.vertList = vertList
        
    def setVertsEdgesFaces(self,vertList:list,edgeList:list,quadList:list):
        self.verts = vertList
        self.edges = edgeList
        self.quads = quadList
    
    def getPosParam(self,xi,eta,zeta):
        N = [(1-xi)*(1-eta)*(1-zeta)/4]
        N.append((1+xi)*(1-eta)*(1-zeta)/4)
        N.append((1+xi)*(1+eta)*(1-zeta)/4)
        N.append((1-xi)*(1+eta)*(1-zeta)/4)
        N.append((1-xi)*(1-eta)*(1+zeta)/4)
        N.append((1+xi)*(1-eta)*(1+zeta)/4)
        N.append((1+xi)*(1+eta)*(1+zeta)/4)
        N.append((1-xi)*(1+eta)*(1+zeta)/4)
        
        cc = 0
        for i in range(8):
            cc += self.verts[i].getCoordinates()*N[i]
        
        return cc
            
    def getQuadList(self):
        return self.quadList
            
    def plotProcessedVertes(self,ax):
        for i in self.edges:
            i.plotProcessedVertes(ax)

        for i in self.vertList:
            vert = self.fem.ptr.getGeometryData().getVertex(i)
            [x,y,z] = vert.getCoordinates()
            ax.scatter(x,y,z,marker='o',color='r')
                
        for i in self.quadList:
            self.plotQuadLineLine(i,ax)
    
    def process(self):
        for i in self.quads:
            i.process()


        pass