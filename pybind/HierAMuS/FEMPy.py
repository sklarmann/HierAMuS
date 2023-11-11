# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

#from CPPFEMPython import HierAMuSPyFEM
#import cppyy

#from FEMPy.pyPointerCollection import pyPointerCollection
from HierAMuS.macroCommands.macroCommands import macroCommands
from HierAMuS.meshCommands.meshCommands import meshCommands
from HierAMuS.plotCommands.plotCommands import plotCommands
from HierAMuS.geomBuilder.geomBuilder import geomBuilder

#from FEMPy.plot.pyVTKplot import pyVTKplot
import os, os.path
import csv

from HierAMuS.HierAMuSPyWrapper import HierAMuSPyFEM

class FEMPy:
    def __init__(self,filepath : str, filename : str):
        self.Info = HierAMuSPyFEM.InfoData()

        # Energy norm and residual norm
        self.maxEnergy = 0
        self.maxResidual = 0

        self.ptr = HierAMuSPyFEM.PointerCollection()
        #self.ptr = pyPointerCollection()
        self.ptr.renew()
        self.ptr.setInfoData(self.Info)
        self.solution = 0
        #self.pyPlot = pyVTKplot()
        #self.ptr.setExternalVtkPlot(self.pyPlot)

        self.propFunctionList = []

        self.SMESH = 0
        self.smesh = 0

        self.csvData = {}


        ff = filename + ".log"
        self.Info.setInFile(ff)
        self.Info.setOutFile(ff)

        ofilename = os.path.join(filepath,ff)
        self.Info.Log.openLogFile(filepath,ff,True)
        self.Info.setInFile(filename)
        self.Info.setOutFile(filename + ".log")
        fname = filepath.replace("\\","/")
        fname += '/'
        self.Info.setDirectory(fname)
        
        self.ptr.getPlotControlInterface().initialize(self.ptr)

        pass

    def __del__(self):
        self.Info.Log.closeLogFile()

    def setMaxThreads(self,n):
        self.ptr.setMaxThreads(n)

    def reset(self):
        self.ptr.renew()

    def setDebug(self,flag):
        self.ptr.setDebug(True)

    def output(self,printLevel : HierAMuSPyFEM.LogLevel, writeLevel : HierAMuSPyFEM.LogLevel, message : str):
        self.Info.outputString(printLevel,writeLevel,message)


    # def setSolutionState(self):
    #     list = HierAMuSPyFEM.ParameterList()
    #     solver=HierAMuSPyFEM.SolverTypes.TypeEigenSimplicialLDLT

    #     solver = HierAMuSPyFEM.StaticSolutionStateHomogenization(self.ptr,list)
    #     self.ptr.setSolutionState(solver)
    #     solver.setSolver(HierAMuSPyFEM.SolverTypes.TypeEigenSimplicialLDLT)

    def setStaticSolutionState(self):
        list = HierAMuSPyFEM.ParameterList()
        solution = HierAMuSPyFEM.StaticSolutionState(list)
        solution.setSolver(HierAMuSPyFEM.SolverTypes.TypeEigenSparseLU)
        self.ptr.setSolutionState(solution)
        
    def setStaticHomogenizationSolutionState(self):
        list = HierAMuSPyFEM.ParameterList()
        solution = HierAMuSPyFEM.StaticSolutionStateHomogenization(list)
        solution.setSolver(HierAMuSPyFEM.SolverTypes.TypeEigenSparseLU)
        self.ptr.setSolutionState(solution)
        
    def getSolutionState(self):
        return self.ptr.getSolutionState()
        

    def setSolver(self,solverNumber):
        """! \brief Set the solver to use.
        \param solverNumber The number of the solver to use.
            1: Eigen SimplicialLLT
            2: Eigen SimplicialLDLT
            3: Eigen SparseLU
            4: Eigen PardisoLDLT
            5: Eigen PardisoLLT
            6: Eigen PardisoLU
        """
        switcher = {
            1: HierAMuSPyFEM.SolverTypes.TypeEigenSimplicialLLT,
            2: HierAMuSPyFEM.SolverTypes.TypeEigenSimplicialLDLT,
            3: HierAMuSPyFEM.SolverTypes.TypeEigenSparseLU,
            4: HierAMuSPyFEM.SolverTypes.TypeEigenPardisoLDLT,
            5: HierAMuSPyFEM.SolverTypes.TypeEigenPardisoLLT,
            6: HierAMuSPyFEM.SolverTypes.TypeEigenPardisoLU
        }

        solver = switcher.get(solverNumber, "Invalid solver")
        self.ptr.getSolutionState().setSolver(solver)


    def getInfoData(self):
        #return self.Info
        pass

    def getMacroCommands(self):
        return macroCommands(self)


    def getMeshCommands(self):
        return meshCommands(self)


    def getPlotCommands(self):
        return plotCommands(self)
        #pass


    def NoLog(self):
        return HierAMuSPyFEM.LogLevel.NoLog

    def BasicLog(self):
        return HierAMuSPyFEM.LogLevel.BasicLog

    def FullLog(self):
        return HierAMuSPyFEM.LogLevel.FullLog

    def outputLine(self,printLevel, writeLevel, output : list):
        pp = ""
        for i in output:
            if type(i) == str:
                pp += i
            else:
                pp += str(i)

        self.Info.Log.outputString(printLevel,writeLevel,pp)

    def basicOutput(self,output : list):
        pp = ""
        for i in output:
            if type(i) == str:
                pp += i
            else:
                pp += str(i)

        self.Info.Log.outputStringBasic(pp)
        
    def debugOutput(self,output : list):
        pp = ""
        for i in output:
            if type(i) == str:
                pp += i
            else:
                pp += str(i)

        self.Info.Log.outputStringDebug(pp)

    def initSalome(self,_SMESH,_smesh):
        self.SMESH = _SMESH
        self.smesh = _smesh

    def addGeomFromSalome(self,mesh):
        b=self.getMeshCommands().getFromSalome()
        b.geomFromSalome(mesh)
        pass

    def getTotalDof(self):
        return self.ptr.getEquationHandler().getTotalDofs()

    def getNumberOfEquations(self):
        return self.ptr.getEquationHandler().getNumberOfEquations()

    def getNumberOfElements(self):
        return self.ptr.getElementList().getNumberOfElements()

    def getFilePath(self):
        return self.Info.getDir()

    def addCSVData(self,name,value):
        self.csvData[name]=value

    def writeCSVData(self,filename,overwrite=False):
        dirName = self.Info.getDir()
        fileName = os.path.join(dirName,filename)

        header = []
        data = []
        for i,v in self.csvData.items():
            header.append(i)
            data.append(v)

        if os.path.isfile(fileName):
            if overwrite:
                ff = open(fileName,'w')
            else:
                ff = open(fileName,'a')

        else:
            ff = open(fileName,'w')
            overwrite=True

        csvWriter = csv.writer(ff)
        if overwrite:
            csvWriter.writerow(header)

        csvWriter.writerow(data)

    def getGeomBuilder(self):
        return geomBuilder(self)



