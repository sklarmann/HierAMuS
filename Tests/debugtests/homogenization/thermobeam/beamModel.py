# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import os,sys
import HierAMuS
import gmsh


L=20
nx=500

gmsh.initialize()
gmsh.model.add("test")

currPath = os.path.dirname(os.path.realpath(__file__))
geoFile = os.path.join(currPath,"RVE.geo")


gmsh.model.occ.addPoint(0,0,0)
gmsh.model.occ.addPoint(L,0,0)

gmsh.model.occ.addLine(1,2)

gmsh.model.occ.synchronize()

gmsh.model.mesh.setTransfiniteCurve(1,nx+1)

gmsh.model.mesh.generate(1)

#gmsh.fltk.run()


fesys = HierAMuS.FEMPy(currPath,"beamModel")
fesys.setStaticSolutionState()
fesys.setSolver(6)
fesys.getMacroCommands().setLogLevel(fesys.FullLog(),fesys.FullLog())

mesh = fesys.getMeshCommands()
macro = fesys.getMacroCommands()
gm = mesh.getFromGMESH()
geo = mesh.getGeometryCommands()


gm.addGeomFromGmsh(gmsh)
geo.checkGeometry()

gm.addLineElements(gmsh,1,1)

mesh.getElementFormulations().addEL103_Timoshenko3D(num=1,E=100,G=50,A=1,Ix=1,Iy=1/12,Iz=1/12,ky=5/6,kz=5/6,geo=0,thermal=1,fe2=0,alpha=1,kappa=1,meshidDisp=1,meshidRot=2,meshidTemp=3,dispOrder=3,rotOrder=3,tempOrder=3)
mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(number=1,E=100,nu=0.3)
mesh.addMaterial(matNum=1,matFormNum=1,elemFormNum=1)

mesh.setDegreesOfFreedom()

mesh.getBoundaryConditions().BCVertex(number=[1],meshId=1,dofs=[1,1,1])
mesh.getBoundaryConditions().BCVertex(number=1,meshId=2,dofs=[1,1,1])
mesh.getBoundaryConditions().BCVertex(number=1,meshId=3,dofs=[1,1,1])
mesh.getBoundaryConditions().BCVertex(number=2,meshId=3,dofs=[1,1,1])


mesh.getBoundaryConditions().BCVertex(number=2,meshId=1,dofs=[1,0,0])
mesh.getBoundaryConditions().BCVertex(number=2,meshId=2,dofs=[1,1,1])

#mesh.getBoundaryConditions().LoadVertex(number=2,meshId=3,load=[1,0,0],propnum=1)
mesh.getBoundaryConditions().PSVertex(number=2,meshId=3,load=[0,0,-2.11],propnum=1)

macro.sparseSetUp()

macro.setPropFunction(1)
macro.setDt(1)

macro.timeincr()

macro.newton(refResidual=1e-9)

fesys.getPlotCommands().toFile()

#macro.printInfo()

E=100
I=1

F=1
print("ref",F*L**3/E/I/3+F*L/50)

print(macro.getSolution(geo.vertexType(),2,1))