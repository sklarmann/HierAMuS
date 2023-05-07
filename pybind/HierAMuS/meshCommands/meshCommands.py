# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from HierAMuS.HierAMuSPyWrapper import HierAMuSPyFEM
from HierAMuS.meshCommands.geometryObjects import geometryObjects
from HierAMuS.meshCommands.finiteElements import finiteElements
from HierAMuS.meshCommands.elementFormulation import elementFormulation
from HierAMuS.meshCommands.materialFormulation import materialFormulation
from HierAMuS.meshCommands.boundaryConditions import boundaryConditions
from HierAMuS.meshCommands.fromSalome import fromSalome
from HierAMuS.meshCommands.fromGMSH import fromGMSh
from HierAMuS.meshCommands.constraintCommands import constraintCommands


from timeit import default_timer as timer

class meshCommands:
    def __init__(self,program):
        self.program = program
        pass

    def getGeometryCommands(self):
        return geometryObjects(self)

    def getElementCommands(self):
        return finiteElements(self)

    def getElementFormulations(self):
        return elementFormulation(self)

    def getMaterialFormulations(self):
        return materialFormulation(self)


    def getBoundaryConditions(self):
        return boundaryConditions(self)

    def getFromSalome(self):
        return fromSalome(self,self.program.SMESH,self.program.smesh)

    def getFromGMESH(self):
        return fromGMSh(self)

    def getConstraintCommands(self):
        return constraintCommands(self)

    def geometryTypes(self):
        return HierAMuSPyFEM.Geometry.GeometryTypes()

    def addMaterial(self,matNum,matFormNum,elemFormNum):
        mat = self.program.ptr.getMaterialList().getMaterial(matNum)
        mat.setElementForumaltion(elemFormNum)
        mat.setMaterialFormulation(matFormNum)

    def setDegreesOfFreedom(self):
        start = timer()
        elemlist = self.program.ptr.getElementList()
        
        elemlist.setDegreesOfFreedom(self.program.ptr)
        end = timer()
        self.program.ptr.getEquationHandler().print(self.program.ptr)
        self.program.basicOutput(["\n   Distributing degrees of freedom took: ", end-start, " seconds"])
        



