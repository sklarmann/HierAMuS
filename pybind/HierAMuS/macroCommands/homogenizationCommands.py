# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from HierAMuS.HierAMuSPyWrapper import HierAMuSPyFEM

from timeit import default_timer as timer

#from HierAMuS.FEMPy import FEMPy
#import HierAMuS

class homogenizationCommands:
    def __init__(self,program):
        self.program = program
        pass
    
    def setHomogenizationSolid2D(self, meshIdDisp, dispOrder, bctype):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("meshIdDisp",meshIdDisp)
        paramList.add("bctype",bctype)
        paramList.add("dispOrder",dispOrder)
        
        self.program.ptr.getSolutionState().initHomogenization(self.program.ptr,1,paramList)
    
    def setHomogenizationBeam(self, meshIdDisp, dispOrder, bctype):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("meshIdDisp",meshIdDisp)
        paramList.add("bctype",bctype)
        paramList.add("dispOrder",dispOrder)
        
        self.program.ptr.getSolutionState().initHomogenization(self.program.ptr,2,paramList)
   
    def setHomogenizationSolid3D(self, meshIdDisp, dispOrder, bctype):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("meshIdDisp",meshIdDisp)
        paramList.add("bctype",bctype)
        paramList.add("dispOrder",dispOrder)
        
        self.program.ptr.getSolutionState().initHomogenization(self.program.ptr,3,paramList)

    def setHomogenizationShell(self, meshIdDisp, dispOrder, bctype):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("meshIdDisp",meshIdDisp)
        paramList.add("bctype",bctype)
        paramList.add("dispOrder",dispOrder)
        
        self.program.ptr.getSolutionState().initHomogenization(self.program.ptr,4,paramList)

    def Homogenization3DThermoMechBeam(self, meshIdDisp, dispOrder, meshIdTemp, tempOrder, bctype):
        paramList = HierAMuSPyFEM.ParameterList()
        
        paramList.add("meshIdDisp",meshIdDisp)
        paramList.add("dispOrder",dispOrder)
        
        paramList.add("meshIdTemp",meshIdTemp)
        paramList.add("tempOrder",tempOrder)
        
        paramList.add("bctype",bctype)
        
        self.program.ptr.getSolutionState().initHomogenization(self.program.ptr,20,paramList)
        
        
    def computeAMatrix(self):
        self.program.ptr.getSolutionState().computeAMatrix(self.program.ptr)
        
    def homogenize(self):
        self.program.ptr.getSolutionState().homogenize(self.program.ptr)

    def setStrains(self,strains : list):
        self.program.ptr.getSolutionState().setStrains(self.program.ptr,strains)
        
    def getCMatrix(self):
        return self.program.ptr.getSolutionState().getCMatrix()

    def getStresses(self):
        return self.program.ptr.getSolutionState().getStresses()
