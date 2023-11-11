# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from HierAMuS.HierAMuSPyWrapper import HierAMuSPyFEM
class boundaryConditions:
    def __init__(self,meshcmd):
        self.mesh = meshcmd
        pass

    def shapeTypes(self):
        return HierAMuSPyFEM.Geometry.ShapeFunctionTypes()
        

    def singleBC(self,eltype,number,meshId,dofs,shapeOrder,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1,set=True):
        self.mesh.program.debugOutput(["WARNING: In meshCommands.boundaryConditions function singleBC is deprecated and replaced by BC"])
        if type(number) == list:
            for i in number:
                elem = self.mesh.program.ptr.getGeometryData().getGeometryElement(eltype,i)
                elem.setBoundaryCondition(meshId,shapeOrder,shapeType,dofs,set)
        else:
            elem = self.mesh.program.ptr.getGeometryData().getGeometryElement(eltype,number)
            elem.setBoundaryCondition(meshId,shapeOrder,shapeType,dofs,set)

    def BC(self,eltype,number,meshId,dofs,shapeOrder,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1,set=True):
        
        if type(number) == list:
            for i in number:
                elem = self.mesh.program.ptr.getGeometryData().getGeometryElement(eltype,i)
                elem.setBoundaryCondition(meshId,shapeOrder,shapeType,dofs,set)
        else:
            elem = self.mesh.program.ptr.getGeometryData().getGeometryElement(eltype,number)
            elem.setBoundaryCondition(meshId,shapeOrder,shapeType,dofs,set)
    
    def BCVertex(self,number,meshId,dofs,set=True):
        ll=[]
        if type(number) == list:
            ll=number
        else:
            ll.append(number)
            
        for i in ll:
            self.BC(self.mesh.getGeometryCommands().vertexType(),i,meshId,dofs,1,HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1,set)
            
    def BCEdge(self,number,meshId,dofs,shapeOrder,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1,set=True):
        self.BC(self.mesh.getGeometryCommands().edgeType(),number,meshId,dofs,shapeOrder,shapeType,set)
        
    def BCFace(self,number,meshId,dofs,shapeOrder,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1,set=True):
        self.BC(self.mesh.getGeometryCommands().faceType(),number,meshId,dofs,shapeOrder,shapeType,set)
        
    def BCVolume(self,number,meshId,dofs,shapeOrder,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1,set=True):
        self.BC(self.mesh.getGeometryCommands().volumeType(),number,meshId,dofs,shapeOrder,shapeType,set)
            
                    
            

    def singleLoad(self,eltype,number,meshId,load,propnum,add=True,localLoad=False,shapeorder=1,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1,direction=[1,0,0]):
        self.mesh.program.debugOutput(["WARNING: In meshCommands.boundaryConditions function singleLoad is deprecated and replaced by load"])
        if type(number) == list:
            for i in number:
                elem = self.mesh.program.ptr.getGeometryData().getGeometryElement(eltype,i)
                elem.setLoad(self.mesh.program.ptr.getLoadList(),meshId,shapeType,shapeorder,load,propnum,direction,localLoad,add)
            
        else:
            elem = self.mesh.program.ptr.getGeometryData().getGeometryElement(eltype,number)
            elem.setLoad(self.mesh.program.ptr.getLoadList(),meshId,shapeType,shapeorder,load,propnum,direction,localLoad,add)
            
    
    def load(self,eltype,number,meshId,load,propnum,add=True,localLoad=False,shapeorder=1,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1,direction=[1,0,0]):
        
        if type(number) == list:
            for i in number:
                elem = self.mesh.program.ptr.getGeometryData().getGeometryElement(eltype,i)
                elem.setLoad(self.mesh.program.ptr.getLoadList(),meshId,shapeType,shapeorder,load,propnum,direction,localLoad,add)
            
        else:
            elem = self.mesh.program.ptr.getGeometryData().getGeometryElement(eltype,number)
            elem.setLoad(self.mesh.program.ptr.getLoadList(),meshId,shapeType,shapeorder,load,propnum,direction,localLoad,add)
    
    def Load(self,eltype,number,meshId,load,propnum,add=True,localLoad=False,shapeorder=1,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1,direction=[1,0,0]):
        ll = []
        if type(number) == list:
            ll = number
        else:
            ll.append(number)
            
            
        for i in ll:
            elem = self.mesh.program.ptr.getGeometryData().getGeometryElement(eltype,i)
            elem.setLoad(self.mesh.program.ptr.getLoadList(),meshId,shapeType,shapeorder,load,propnum,direction,localLoad,add)
            
    def LoadVertex(self,number,meshId,load,propnum,add=True):
        self.Load(self.mesh.getGeometryCommands().vertexType(),number,meshId,load,propnum,add)
        
    def LoadEdge(self,number,meshId,load,propnum,add=True,shapeorder=1,localload=False,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1,direction=[1,0,0]):
        self.Load(self.mesh.getGeometryCommands().edgeType(),number,meshId,load,propnum,add,shapeorder=shapeorder,shapeType=shapeType,direction=direction)
        
    def LoadFace(self,number,meshId,load,propnum,add=True,shapeorder=1,localload=False,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1,direction=[1,0,0]):
        self.Load(self.mesh.getGeometryCommands().faceType(),number,meshId,load,propnum,add,shapeorder=shapeorder,shapeType=shapeType,direction=direction)
        
    def LoadVolume(self,number,meshId,load,propnum,add=True,shapeorder=1,localload=False,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1,direction=[1,0,0]):
        self.Load(self.mesh.getGeometryCommands().volumeType(),number,meshId,load,propnum,add,shapeorder=shapeorder,shapeType=shapeType,direction=direction)
    
    
    def prescribedSolution(self,eltype,number,meshId,load,propnum,add=True,localLoad=False,shapeorder=1,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1,direction=[1,0,0]):
        
        if type(number) == list:
            for i in number:
                elem = self.mesh.program.ptr.getGeometryData().getGeometryElement(eltype,i)
                elem.setLoad(self.mesh.program.ptr,meshId,shapeType,shapeorder,load,propnum,direction,localLoad,add)
            
        else:
            elem = self.mesh.program.ptr.getGeometryData().getGeometryElement(eltype,number)
            elem.setLoad(self.mesh.program.ptr,meshId,shapeType,shapeorder,load,propnum,direction,localLoad,add)
            
            
    def PS(self,eltype,number,meshId,load,propnum,add=True,localLoad=False,shapeorder=1,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1,direction=[1,0,0]):
        
        ll=[]
        if type(number) == list:
            ll = number
        else:
            ll.append(number)
        
        for i in ll:
            elem = self.mesh.program.ptr.getGeometryData().getGeometryElement(eltype,i)
            elem.setPrescribedSolution(self.mesh.program.ptr,meshId,shapeType,shapeorder,load,propnum,direction,localLoad,add)
            
    
    def PSVertex(self,number,meshId,load,propnum,add=True):
        self.PS(self.mesh.getGeometryCommands().vertexType(),number,meshId,load,propnum,add)
        
    def PSFaceH1(self,number,meshId,load,propnum,add=True,localLoad=False,shapeorder=1,direction=[1,0,0]):
        
        self.PS(self.mesh.getGeometryCommands().faceType(),number,meshId,load,propnum,add=add,localLoad=localLoad,shapeorder=shapeorder,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1,direction=[1,0,0])






