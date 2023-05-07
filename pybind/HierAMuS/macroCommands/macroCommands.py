# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from HierAMuS.HierAMuSPyWrapper import HierAMuSPyFEM
from HierAMuS.macroCommands.homogenizationCommands import homogenizationCommands

from timeit import default_timer as timer
from string import Template

#from HierAMuS.FEMPy import FEMPy
#import HierAMuS


def NameTimeString(name,timeval):
    return "{:40}{:8.4f} sec".format(name,timeval)

def NameFloatString(name,timeval):
    return "{:40}{:12.6e}".format(name,timeval)

class macroCommands:
    def __init__(self,program):
        self.program = program
        pass
    
    def getHomogenizationCommands(self):
        return homogenizationCommands(self.program)


    def setLogLevel(self,printNum, writeNum):
        self.program.Info.Log.setLogLevel(printNum,writeNum)

    def printInfo(self):
        self.program.ptr.getGeometryData().print(self.program.ptr)
        self.program.ptr.getLoadList().print(self.program.ptr)
        props = self.program.ptr.getPropLoads()
        props.print(self.program.ptr)

    def sparseSetUp(self):
        # elemList = self.program.ptr.getElementList()
        # numElems = elemList.getNumberOfElements()
        # self.program.ptr.getEquationHandler().updateEquations()
        # self.program.ptr.getEquationHandler().initSolutionState()

        # for i in range(numElems):
        #     elemList.getElement(i).setUpSparseMatrix(self.program.ptr)

        start = timer()
        self.program.ptr.getSolutionState().setSparseMatrix(self.program.ptr)
        end = timer()
        self.program.ptr.getEquationHandler().print(self.program.ptr)
        
        self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),NameTimeString("Setting up Sparse Matrix took:",end-start))
        
    def computeNorm(self,param):
        normType = HierAMuSPyFEM.FiniteElement.NormTypes.L2Stresses
        elemList = self.program.ptr.getElementList()
        numElems = elemList.getNumberOfElements()

        res = 0
        for i in range(numElems):
            res += elemList.getElement(i).computeNorm(self.program.ptr,normType)

        return res

    # def setTimeName(self,timename="t",dtName="dt"):
    #     tname = cppyy.gbl.std.string(timename)
    #     dtname = cppyy.gbl.std.string(dtName)

    #     self.program.ptr.getSolutionState().getProps().setNames(self.program.ptr,tname,dtname)

    #     #print(self.program.ptr.getPropLoads())

    def setPropFunction(self,number,function=lambda t : t,tmin=0,tmax=100):
        #self.program.propFunctionList.append(function)
        a = HierAMuSPyFEM.Function()
        a.setLambda(function)
        self.program.ptr.getSolutionState().getProps().addFunction(self.program.ptr,number,a,tmin,tmax)

    def setDt(self,dt):
        self.program.ptr.getSolutionState().getProps().set_dt(dt)

    def timeincr(self):
        self.program.ptr.getSolutionState().nextSolutionStep()
        self.program.ptr.getSolutionState().updateRVEHistory(self.program.ptr)
        t=self.program.ptr.getPropLoads().getTime()
        self.program.ptr.getPlotControlInterface().timeUpdate(t)
        
        self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),"Time Update to {:>12}{:12.6e}".format("t=",t))


    def assembleSolve(self):
        start = timer()
        soltstate = self.program.ptr.getSolutionState()
        soltstate.assembleSystem(self.program.ptr)
        end = timer()
        
        self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),NameTimeString("Assembly of equation system took:",end-start))

        soltstate.factorize()
        start = timer()
        self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),NameTimeString("Factorization of equation system took:",start-end))

        soltstate.solve(self.program.ptr)
        soltstate.updateSolution(self.program.ptr)
        end = timer()
        self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),NameTimeString("Solution of equation system took:",end-start))
        energy0 = soltstate.energyNorm()
        residual0 = soltstate.residual()
        self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),NameFloatString("Initial Residual norm was:",residual0))

        #self.program.ptr.solution.solve()

    def newton(self,maxIteration=10,refResidual=1e-13,damped=False):
        iteration = 1
        self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),["------------------------------------------------------------------------"])
        self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),"Starting newton iterations")

        solutionState = self.program.ptr.getSolutionState()
        start = timer()
        self.assembleSolve()
        end = timer()
        st = "{:15}{:6d}/{:6d}  {:9}:{:8.4f} sec".format("Iteration",iteration,maxIteration,"took",end-start)
        self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),st)
        energy0 = solutionState.energyNorm()
        if energy0 == 0:
            energy0 = 1
        residual0 = solutionState.residual()
        if residual0 == 0:
            residual0 = 1

        energyNom = 1
        residualNorm = 1

        
        self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),NameFloatString("Current residual:",residual0))
        self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),NameFloatString("Current energy::",energy0))
        self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),NameFloatString("Residual norm::",residualNorm))
        self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),NameFloatString("Energy Norm:",energyNom))
        self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),"------------------------------------------------------------------------")

        while iteration < maxIteration and (residualNorm > refResidual or energyNom > refResidual):
            iteration+=1
            start = timer()
            self.assembleSolve()
            end = timer()
            energy = solutionState.energyNorm()
            residual = solutionState.residual()
            residualNorm = residual/residual0
            energyNom = energy/energy0

            
            st = "{:15}{:6d}/{:6d}  {:9}:{:8.4f} sec".format("Iteration",iteration,maxIteration,"took",end-start)
            self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),st)
            self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),NameFloatString("Current residual:",residual))
            self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),NameFloatString("Current energy::",energy))
            self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),NameFloatString("Residual norm::",residualNorm))
            self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),NameFloatString("Energy Norm:",energyNom))
            self.program.outputLine(self.program.BasicLog(),self.program.BasicLog(),"------------------------------------------------------------------------")

    def getSolution(self,geomType,geomNumber,meshId):
        geoData = self.program.ptr.getGeometryData()
        solState = self.program.ptr.getSolutionState()
        geo = geoData.getGeometryElement(geomType,geomNumber)
        nodes=geo.getNodes(self.program.ptr,meshId)

        sol = []
        for nn in nodes:
            dofs = nn.getDegreesOfFreedom(self.program.ptr)
            tsol = solState.getSolution(dofs)
            for ss in tsol:
                sol.append(ss)

        return sol

    def computeEigenValues(self,numOfEV,numOfSearchEV=0,max=False, convTolerance=1e-10, shift = 1e-10):
        solState = self.program.ptr.getSolutionState()
        solState.computeEigenValues(self.program.ptr,numOfEV,numOfSearchEV,max,convTolerance,shift)

    def printSpMat(self):
        solState = self.program.ptr.getSolutionState()
        solState.printSpMat()
        
    def setStrains(self,strains):
        sol = self.program.ptr.getSolutionState()
        sol.setStrains(strains)
        
    def solutionStateToFile(self,filename):
        self.program.ptr.solutionStateToFile(filename)
        
    def solutionStateFromFile(self,filename):
        self.program.ptr.solutionStateFromFile(filename)