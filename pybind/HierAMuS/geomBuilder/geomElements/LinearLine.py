# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from HierAMuS.geomBuilder.geomElements.Vertex import Vertex
import numpy as np

class Line:
    def __init__(self,fesys):
        self.geoId = -1
        self.fem = fesys
        self.ndiv = 0
        self.vertList = []
        self.EdgeList = []
        pass

    def setNumElements(self,ndiv):
        if self.ndiv == 0:
            self.ndiv = ndiv
        elif self.ndiv != ndiv:
            print('Error, tried to change number of divisions for Line ', self.geoId, ' changed')
            

    def setElementSize(self,esize):
        dx = (self.v2.getCoordinates()-self.v1.getCoordinates())
        L = np.linalg.norm(dx)
        div = round(L/esize)
        if self.ndiv==0:
            self.ndiv = int(L/esize)
        elif self.ndiv != div:
            print('Error, tried to change number of divisions for Line ', self.geoId, ' changed')
            
    def getLength(self):
        dx = (self.v2.getCoordinates()-self.v1.getCoordinates())
        L = np.linalg.norm(dx)
        return L
        

    def getVertList(self,orient):
        if orient == 1:
            return self.vertList
        else:
            t = self.vertList.copy()
            t.reverse()
            return t

    def getEdgeList(self,orient):
        if orient == 1:
            return self.EdgeList
        else:
            t = self.EdgeList.copy()
            t.reverse()
            return t

    def setStartEnd(self,xstart : Vertex, xend : Vertex):
        self.v1 = xstart
        self.v2 = xend

    def getPosParam(self,xi):
        c1 = self.v1.getCoordinates()
        c2 = self.v2.getCoordinates()
        
        ccoor = c1*(1-xi)/2+c2*(1+xi)/2
        return ccoor
    
    def getOrientation(self,vert:Vertex):
        if vert.geoId == self.v1.geoId:
            return 1
        else:
            return -1

    def process(self):
        dx = (self.v2.getCoordinates()-self.v1.getCoordinates())/self.ndiv

        self.vertList.append(self.v1.num)

        for i in range(self.ndiv-1):
            initCoor = self.v1.getCoordinates()
            x = initCoor[0] + (i+1)*dx[0]
            y = initCoor[1] + (i+1)*dx[1]
            z = initCoor[2] + (i+1)*dx[2]
            #print(x,y,z)
            self.vertList.append(self.fem.getMeshCommands().getGeometryCommands().requestVertex(x, y, z))

        self.vertList.append(self.v2.num)

        for i in range(self.ndiv):
            vl = self.vertList[i:i+2]
            self.EdgeList.append(self.fem.getMeshCommands().getGeometryCommands().requestLinearEdge(vl))


        pass

    def plotSelfInitial(self,ax):
        coora = self.v1.getCoordinates()
        coorb = self.v2.getCoordinates()
        x = [coora[0],coorb[0]]
        y = [coora[1],coorb[1]]
        z = [coora[2],coorb[2]]
        ax.plot(x,y,z,color='b')

    def plotProcessedVertes(self,ax):
        for i in self.vertList:
            vert = self.fem.ptr.getGeometryData().getVertex(i)
            coor = vert.getCoordinates()

            ax.scatter(coor[0],coor[1],coor[2],marker='o',color='r')
            