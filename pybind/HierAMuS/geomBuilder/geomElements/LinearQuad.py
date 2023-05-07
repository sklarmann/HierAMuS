# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from HierAMuS.geomBuilder.geomElements.Vertex import Vertex
from HierAMuS.geomBuilder.geomElements.LinearLine import Line

class Quad:
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


    def setElementSize(self,esx,esy):
        Lx1 = self.edges[0].getLength()
        Lx2 = self.edges[2].getLength()
        Ly1 = self.edges[1].getLength()
        Ly2 = self.edges[3].getLength()
        
        Lax = (Lx1+Lx2)/2
        Lay = (Ly1+Ly2)/2
        
        ndivx = round(Lax/esx)
        ndivy = round(Lay/esy)
        
        self.setDivision(ndivx,ndivy)


    def setDivision(self,nx,ny):
        self.nx = nx
        self.ny = ny

        self.edges[0].setNumElements(nx)
        self.edges[1].setNumElements(ny)
        self.edges[2].setNumElements(nx)
        self.edges[3].setNumElements(ny)
        
    def setVertsEdges(self,vertList:list,edgeList:list):
        self.verts = vertList
        self.edges = edgeList
    
    def getPosParam(self,xi,eta):
        Nxi1 = (1-xi)/2
        Nxi2 = (1+xi)/2
        
        e2 = self.edges[1]
        e4 = self.edges[3]
        
        c2 = e2.getPosParam(eta*e2.getOrientation(self.verts[1]))
        c4 = e4.getPosParam(eta*e4.getOrientation(self.verts[0]))
        
        
        
        
        return Nxi2*c2+Nxi1*c4
            
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
        dxi = 2/self.nx
        deta = 2/self.ny

        for i in range(self.ny-1):
            for j in range(self.nx-1):
                xi = (j+1)*dxi-1
                eta = (i+1)*deta-1
                [x,y,z] = self.getPosParam(xi,eta)
                self.vertList.append(self.fem.getMeshCommands().getGeometryCommands().requestVertex(x, y, z))

        # lower boundary
        edge1 = self.edges[0]
        edge2 = self.edges[1]
        edge3 = self.edges[2]
        edge4 = self.edges[3]
        vertex1 = self.verts[0]
        vertex2 = self.verts[1]
        vertex3 = self.verts[2]
        vertex4 = self.verts[3]

        e1EdgeList = edge1.getEdgeList(edge1.getOrientation(vertex1))
        e2EdgeList = edge2.getEdgeList(edge2.getOrientation(vertex2))
        e3EdgeList = edge3.getEdgeList(edge3.getOrientation(vertex4))
        e4EdgeList = edge4.getEdgeList(edge4.getOrientation(vertex1))


        e1VertList = edge1.getVertList(edge1.getOrientation(vertex1))
        e2VertList = edge2.getVertList(edge2.getOrientation(vertex2))
        e3VertList = edge3.getVertList(edge3.getOrientation(vertex4))
        e4VertList = edge4.getVertList(edge4.getOrientation(vertex1))


        # setup verts list
        verts = e1VertList
        pos = 0
        for j in range(self.ny-1):
            verts += [e4VertList[j+1]]
            for i in range(self.nx-1):
                verts += [self.vertList[pos]]
                pos+=1
            verts += [e2VertList[j+1]]
        verts += e3VertList

        # setup edges list
        hedges = e1EdgeList
        vedges = []
        posl = 0
        posu = self.nx+1
        for j in range(self.ny-1):
            vedges.append(e4EdgeList[j])
            [v1,v2,v3,v4] = [verts[posl],verts[posl+1],verts[posu+1],verts[posu]]
            hedges.append(self.fem.getMeshCommands().getGeometryCommands().requestLinearEdge([v4,v3]))
            for i in range(self.nx - 1):
                posl+=1
                posu+=1
                [v1,v2,v3,v4] = [verts[posl],verts[posl+1],verts[posu+1],verts[posu]]
                vedges.append(self.fem.getMeshCommands().getGeometryCommands().requestLinearEdge([v1,v4]))
                hedges.append(self.fem.getMeshCommands().getGeometryCommands().requestLinearEdge([v4,v3]))

            
            posl+=1
            posu+=1
            [v1,v2,v3,v4] = [verts[posl],verts[posl+1],verts[posu+1],verts[posu]]
            vedges.append(e2EdgeList[j])
            posl+=1
            posu+=1

        hedges += e3EdgeList
        vedges.append(e4EdgeList[self.ny-1])
        for i in range(self.nx - 1):
            posl+=1
            posu+=1
            [v1,v2,v3,v4] = [verts[posl],verts[posl+1],verts[posu+1],verts[posu]]
            vedges.append(self.fem.getMeshCommands().getGeometryCommands().requestLinearEdge([v1,v4]))
        vedges.append(e2EdgeList[self.ny-1])

        posl = 0
        posu = self.nx+1
        posle = 0
        posue = self.nx
        for i in range(self.ny):
            for j in range(self.nx):
                vnums = [verts[posl],verts[posl+1],verts[posu+1],verts[posu]]
                enums = [hedges[posle],vedges[posl+1],hedges[posue],vedges[posl]]
                self.quadList.append(self.fem.getMeshCommands().getGeometryCommands().requestLinearQuadrilateralFace(vnums, enums))
                posl += 1
                posu += 1
                posle+= 1
                posue+= 1

            posl += 1
            posu += 1




        pass