# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

class Vertex:
    def __init__(self,fesys):
        self.geoId = -1
        self.fem = fesys
        self.num = -1
        pass

    def setCoordinates(self,x,y,z):
        if self.num == -1:
            self.num = self.fem.getMeshCommands().getGeometryCommands().requestVertex(x, y, z)
        else:
            vert = self.fem.ptr.getGeometryData().getVertex(self.num)
            vert.setCoordinates(x,y,z)

    def getCoordinates(self):
        vert = self.fem.ptr.getGeometryData().getVertex(self.num)
        return vert.getCoordinates()

    def plotSelf(self,ax):
        coor = self.getCoordinates()
        ax.scatter(coor[0],coor[1],coor[2],marker='o',color='r')
        ax.text(coor[0],coor[1],coor[2],str(self.num))