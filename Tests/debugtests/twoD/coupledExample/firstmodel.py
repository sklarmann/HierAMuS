# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import HierAMuS
import sys, os
import gmsh


path = os.path.dirname(__file__)

# Constructing the RVE
#geom
L=2
h=2
#mesh
nx=10
ny=10
#solutionsettings
meshid=1 
disporder=2

gmsh.initialize()
gmsh.option.setNumber('General.Terminal', 0)
gmsh.model.add("RVE")

p1=gmsh.model.occ.addPoint(-L/2,-h/2,0)
p2=gmsh.model.occ.addPoint(L/2,-h/2,0)
p3=gmsh.model.occ.addPoint(L/2,h/2,0)
p4=gmsh.model.occ.addPoint(-L/2,h/2,0)

l1=gmsh.model.occ.addLine(p1,p2)
l2=gmsh.model.occ.addLine(p2,p3)
l3=gmsh.model.occ.addLine(p4,p3)
l4=gmsh.model.occ.addLine(p1,p4)


cl = gmsh.model.occ.addCurveLoop([l1,l2,l3,l4])
f1 = gmsh.model.occ.addPlaneSurface([cl])

gmsh.model.occ.synchronize()

gmsh.model.mesh.setTransfiniteCurve(l1,nx+1)
gmsh.model.mesh.setTransfiniteCurve(l3,nx+1)
gmsh.model.mesh.setTransfiniteCurve(l2,ny+1)
gmsh.model.mesh.setTransfiniteCurve(l4,ny+1)
gmsh.model.mesh.setTransfiniteSurface(f1)

gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.model.mesh.generate(2)


RVE=HierAMuS.FEMPy(path,"RVE")
RVE.setStaticHomogenizationSolutionState()
RVE.setSolver(4)

gm = RVE.getMeshCommands().getFromGMESH()
geo = RVE.getMeshCommands().getGeometryCommands()
gm.addGeomFromGmsh(gmsh)
geo.checkGeometry()

gm.addFaceElements(gmsh,f1,1)
RVE.getMeshCommands().getElementFormulations().addEL201_2DShell(num=1,meshiddisp=meshid,disporder=disporder,mode=1)
RVE.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(number=1,E=100,nu=0.3,thickness=1,plainstrain=1)
RVE.getMeshCommands().addMaterial(matNum=1,matFormNum=1,elemFormNum=1)

RVE.getMeshCommands().setDegreesOfFreedom()

lline = gm.getEdgeNumbers(gmsh,l4,1)
rline = gm.getEdgeNumbers(gmsh,l2,1)
tline = gm.getEdgeNumbers(gmsh,l3,1)
bline = gm.getEdgeNumbers(gmsh,l1,1)

RVE.getMacroCommands().getHomogenizationCommands().setHomogenizationSolid2D(meshIdDisp=meshid,dispOrder=disporder,bctype=0)

RVE.getMacroCommands().sparseSetUp()
RVE.getMacroCommands().getHomogenizationCommands().computeAMatrix()

RVE.getMacroCommands().getHomogenizationCommands().setStrains([0.0,0.0,0.0])

gmsh.finalize()


# Constructing the Macroscale model
#geom
L=10
h=1
#mesh
nx=2
ny=2
#solutionsettings
meshid=1 
disporder=2

gmsh.initialize()
gmsh.option.setNumber('General.Terminal', 0)
gmsh.model.add("Macro")

p1=gmsh.model.occ.addPoint(-L/2,-h/2,0)
p2=gmsh.model.occ.addPoint(L/2,-h/2,0)
p3=gmsh.model.occ.addPoint(L/2,h/2,0)
p4=gmsh.model.occ.addPoint(-L/2,h/2,0)

l1=gmsh.model.occ.addLine(p1,p2)
l2=gmsh.model.occ.addLine(p2,p3)
l3=gmsh.model.occ.addLine(p4,p3)
l4=gmsh.model.occ.addLine(p1,p4)


cl = gmsh.model.occ.addCurveLoop([l1,l2,l3,l4])
f1 = gmsh.model.occ.addPlaneSurface([cl])

gmsh.model.occ.synchronize()

gmsh.model.mesh.setTransfiniteCurve(l1,nx+1)
gmsh.model.mesh.setTransfiniteCurve(l3,nx+1)
gmsh.model.mesh.setTransfiniteCurve(l2,ny+1)
gmsh.model.mesh.setTransfiniteCurve(l4,ny+1)
gmsh.model.mesh.setTransfiniteSurface(f1)

gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.model.mesh.generate(2)

MACRO=HierAMuS.FEMPy(path,"MACRO")
MACRO.setStaticSolutionState()
MACRO.setSolver(4)
MACRO.getMacroCommands().setLogLevel(MACRO.BasicLog(),MACRO.BasicLog())

gm = MACRO.getMeshCommands().getFromGMESH()
geo = MACRO.getMeshCommands().getGeometryCommands()
gm.addGeomFromGmsh(gmsh)
geo.checkGeometry()

gm.addFaceElements(gmsh,f1,1)
MACRO.getMeshCommands().getElementFormulations().addEL201_2DShell(num=1,meshiddisp=meshid,disporder=disporder,mode=1)
MACRO.getMeshCommands().getMaterialFormulations().addMAS1_Homogenization(number=1,RVE=RVE)
#MACRO.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(number=1,E=100,nu=0.3,thickness=1,plainstrain=1)
MACRO.getMeshCommands().addMaterial(matNum=1,matFormNum=1,elemFormNum=1)

MACRO.getMeshCommands().setDegreesOfFreedom()

bclines = gm.getEdgeNumbers(gmsh,l4,1)
MACRO.getMeshCommands().getBoundaryConditions().BC(eltype=geo.edgeType(),number=bclines,meshId=meshid,dofs=[1,1,1],shapeOrder=disporder)
loadlines = gm.getEdgeNumbers(gmsh,l2,1)
MACRO.getMeshCommands().getBoundaryConditions().load(eltype=geo.edgeType(),number=loadlines,meshId=meshid,load=[0,1,0],propnum=1)

MACRO.getMacroCommands().sparseSetUp()

MACRO.getMacroCommands().setPropFunction(1)
MACRO.getMacroCommands().setDt(1)
MACRO.getMacroCommands().timeincr()

MACRO.getMacroCommands().newton(refResidual=1e-10)

MACRO.getPlotCommands().toFile()

#gmsh.fltk.run()