# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from CPPFEMPython import HierAMuSPyFEM
import cppyy
import paraview.simple
import paraview.vtk
import os
import vtk

class pyVTKplot(HierAMuSPyFEM.vtkPlotInterface):
    
    def __init__(self) -> None:
        self.mainGrid = paraview.vtk.vtkMultiBlockDataSet()
        
        # FieldTypes 
        self.PointDataField = "PointField"
        self.CellDataField = "CellField"
        
        # Paraview Field Names
        self.intName = "int"
        self.precName = "double"
        
        self.pointIdName = "PointId"
        self.pointIdMap = cppyy.gbl.std.map["indexType",cppyy.gbl.std.map["indexType",cppyy.gbl.std.map["indexType","indexType"]]]()
        
        self.FEcellIdField = "Element Number"
        self.FESubCellIdField = "Sub Element Number"
        self.cellIdMap = cppyy.gbl.std.map["indexType",cppyy.gbl.std.map["indexType",cppyy.gbl.std.map["indexType",cppyy.gbl.std.map["indexType","indexType"]]]]()

        self.time = []
        
        super().__init__()

    def timeUpdate(self,time):
        self.time.append(time.value)

    def setFilePathName(self,path,name):
        self.pathName = os.path.join(path,'parvout')
        self.fileName = name
        self.pvdFile = os.path.join(self.pathName,self.fileName + '.pvd')
        

    def writePVD(self):
        if os.path.exists(self.pvdFile):
            os.remove(self.pvdFile)

        f = open(self.pvdFile,'w')
        f.write('<?xml version ="1.0"?>\n')
        f.write('<VTKFile type="Collection" version="0.1">\n')
        f.write('  <Collection>\n')
        for i in range(len(self.time)):
            f.write('    <DataSet timestep="'+ str(self.time[i]) +'" group="" part=""\n')
            filepath = self.fileName + '/' + self.fileName + str(i) + '.vtm'
            f.write('    file="'+filepath+'"/>\n')

        f.write('  </Collection>\n')
        f.write('</VTKFile>')

        

    def getMainBlock(self,main):
        a = self.mainGrid.GetBlock(main)
        if a:
            return a
        else:
            b = paraview.vtk.vtkMultiBlockDataSet()
            self.mainGrid.SetBlock(main,b)
            a = self.mainGrid.GetBlock(main)
            return a

    def getGrid(self,main,part):
        block = self.getMainBlock(main)

        grid = block.GetBlock(part)
        if grid:
            return grid
        else:
            grid = paraview.vtk.vtkUnstructuredGrid()
            block.SetBlock(part,grid)
            pointList = paraview.vtk.vtkPoints()
            grid.SetPoints(pointList)
            return block.GetBlock(part)
        
    def getParaviewPointNumber(self,main,part,femId):
        if self.pointIdMap.find(main) != self.pointIdMap.end():
            if self.pointIdMap[main].find(part) != self.pointIdMap[main].end():
                if self.pointIdMap[main][part].find(femId) != self.pointIdMap[main][part].end():
                    return True,self.pointIdMap[main][part][femId]
                
        parvId = self.pointIdMap[main][part].size()
        self.pointIdMap[main][part][femId] = parvId
        return False,parvId
    
    def getParaviewCellNumber(self,main,part,FEId,FESubId):
        if self.cellIdMap.find(main) != self.cellIdMap.end():
            if self.cellIdMap[main].find(part) != self.cellIdMap[main].end():
                if self.cellIdMap[main][part].find(FEId) != self.cellIdMap[main][part].end():
                    if self.cellIdMap[main][part][FEId].find(FESubId) != self.cellIdMap[main][part][FEId].end():
                        return True,self.pointIdMap[main][part][FEId][FESubId]
                
        grid = self.getGrid(main,part)
        parvId=grid.GetNumberOfCells()
        self.cellIdMap[main][part][FEId][FESubId] = parvId
        return False,parvId
        
    def getField(self,main,part,fieldName,fieldType,numberType,numComp):
        grid = self.getGrid(main,part)
        
        field = 0
        if fieldType == self.PointDataField:
            field = grid.GetPointData().GetArray(fieldName)
        elif fieldType == self.CellDataField:
            field = grid.GetCellData().GetArray(fieldName)
        else:
            print("Error, something wrong in getField!", fieldName)
             
        if not field:
            if numberType == self.intName:
                field = paraview.vtk.vtkIntArray()
            elif numberType == self.precName:
                field = paraview.vtk.vtkDoubleArray()
                
            field.SetName(fieldName)
            field.SetNumberOfComponents(numComp)
            field.SetNumberOfTuples(grid.GetNumberOfPoints())
            if fieldType == self.PointDataField:
                grid.GetPointData().AddArray(field)
            elif fieldType == self.CellDataField:
                grid.GetCellData().AddArray(field)
            #print("Non existiend field added!", fieldName)
        
        
        return field
    
    def writeFile(self):
        if len(self.time) == 0:
            self.time.append(0)
        
        num = len(self.time)
        pathfile=os.path.join(self.pathName,self.fileName + '/' + self.fileName + str(num-1) + '.vtm')
        writer = vtk.vtkXMLMultiBlockDataWriter()
        writer.SetFileName(pathfile)
        writer.SetInputData(self.mainGrid)
        writer.SetDataModeToAscii()
        
        writer.Write()
        self.writePVD()

        
    def addPoint(self,main,part,lId,x,y,z):
        grid = self.getGrid(main,part)
        field = self.getField(main,part,self.pointIdName,self.PointDataField,self.intName,1)
        #grid.GetPoints().InsertPoint(lId,x,y,z)
        present,parvId = self.getParaviewPointNumber(main,part,lId)
        if present:
            # Update Coordinates
            grid.GetPoints().InsertPoint(parvId,x,y,z)
            
        else:
            grid.GetPoints().InsertNextPoint(x,y,z)
            field.InsertComponent(parvId,0,lId)
            
            
        
    def addCell(self,main,part,FEId,FESubId,points,numpts,vtkNum):
        grid = self.getGrid(main,part)
        present,parvId = self.getParaviewCellNumber(main,part,FEId,FESubId)
        
        list = paraview.vtk.vtkIdList()
        for i in range(numpts):
            present, pId = self.getParaviewPointNumber(main,part,points[i])
            if present:
                list.InsertNextId(pId)
            
        feIdField = self.getField(main,part,self.FEcellIdField,self.CellDataField,self.intName,1)
        feSubIdField = self.getField(main,part,self.FESubCellIdField,self.CellDataField,self.intName,1)
        feIdField.InsertComponent(parvId,0,FEId)
        feSubIdField.InsertComponent(parvId,0,FESubId)
        grid.InsertNextCell(vtkNum,list)
        
        
    def setPointData(self,main,part,FEId,data,num_comp,name):
        name = name.c_str()
        field = self.getField(main,part,name,self.PointDataField,self.precName,num_comp)
        present,parvId = self.getParaviewPointNumber(main,part,FEId)
        #print(data)
        for i in range(num_comp):
            field.InsertComponent(parvId,i,data[i])
        
        
    