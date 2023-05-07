# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
from HierAMuS.geomBuilder.geomElements import Vertex, LinearLine, LinearQuad

class geomBuilder:
    def __init__(self,fesys):
        self.fem = fesys
        self.vertList = []
        self.lineList = []
        self.quadList = []
        pass

    def getVertex(self):
        t = Vertex.Vertex(self.fem)
        t.geoId = len(self.vertList)
        self.vertList.append(t)
        return t

    def getEdge(self):
        l = LinearLine.Line(self.fem)
        l.geoId = len(self.lineList)
        self.lineList.append(l)
        return l
    
    def getQuad(self):
        q = LinearQuad.Quad(self.fem)
        q.geoId = len(self.quadList)
        self.quadList.append(q)
        return q

    def process(self):

        for i in self.lineList:
            i.process()
            
        for i in self.quadList:
            i.process()

    def plot(self,fig):
        ax = fig.add_subplot(221,projection='3d')
        ax2 = fig.add_subplot(222,projection='3d')
        ax3 = fig.add_subplot(223,projection='3d')
        ax4 = fig.add_subplot(224,projection='3d')

        ax.set_title('Initial Points')
        ax2.set_title('Initial Edges')
        ax3.set_title('Generated Points')
        ax4.set_title('Generated Points')

        for i in self.vertList:
            i.plotSelf(ax)

        for i in self.lineList:
            i.plotSelfInitial(ax2)

        for i in self.lineList:
            i.plotProcessedVertes(ax3)

        for i in self.quadList:
            i.plotProcessedVertes(ax4)
