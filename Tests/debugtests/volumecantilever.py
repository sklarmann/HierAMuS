# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import gmsh
import sys, os
import HierAMuS
import numpy as np

gmsh.initialize()

gmsh.model.add("t11")

L=10
b=1
h=1

nx=4
ny=4
nz=4

disporder=2

E=100
nu=0.3
fy=10
xh=E/4
e=1

p1 = gmsh.model.occ.addPoint(0, 0, 0)
p2 = gmsh.model.occ.addPoint(L, 0, 0)
p3 = gmsh.model.occ.addPoint(L, b, 0)
p4 = gmsh.model.occ.addPoint(0, b, 0)
p5 = gmsh.model.occ.addPoint(0, 0, h)
p6 = gmsh.model.occ.addPoint(L, 0, h)
p7 = gmsh.model.occ.addPoint(L, b, h)
p8 = gmsh.model.occ.addPoint(0, b, h)

l1=gmsh.model.occ.addLine(p1, p2)
l2=gmsh.model.occ.addLine(p2, p3)
l3=gmsh.model.occ.addLine(p3, p4)
l4=gmsh.model.occ.addLine(p4, p1)

l5=gmsh.model.occ.addLine(p1, p5)
l6=gmsh.model.occ.addLine(p2, p6)
l7=gmsh.model.occ.addLine(p3, p7)
l8=gmsh.model.occ.addLine(p4, p8)

l9 =gmsh.model.occ.addLine(p5, p6)
l10=gmsh.model.occ.addLine(p6, p7)
l11=gmsh.model.occ.addLine(p7, p8)
l12=gmsh.model.occ.addLine(p8, p5)

cl1 = gmsh.model.occ.addCurveLoop([l1,l2,l3,l4])
f1 = gmsh.model.occ.addPlaneSurface([cl1])
cl2 = gmsh.model.occ.addCurveLoop([l1,l6,l9,l5])
f2 = gmsh.model.occ.addPlaneSurface([cl2])
cl3 = gmsh.model.occ.addCurveLoop([l2,l7,l10,l6])
f3 = gmsh.model.occ.addPlaneSurface([cl3])
cl4 = gmsh.model.occ.addCurveLoop([l3,l8,l11,l7])
f4 = gmsh.model.occ.addPlaneSurface([cl4])
cl5 = gmsh.model.occ.addCurveLoop([l4,l5,l12,l8])
f5 = gmsh.model.occ.addPlaneSurface([cl5])
cl6 = gmsh.model.occ.addCurveLoop([l9,l10,l11,l12])
f6 = gmsh.model.occ.addPlaneSurface([cl6])

fl = gmsh.model.occ.addSurfaceLoop([f1,f2,f3,f4,f5,f6])
vol = gmsh.model.occ.addVolume([fl])
#curveLoop = gmsh.model.occ.addCurveLoop([-l1, -l2, l3, -l4, -circ])
#surface = gmsh.model.occ.addPlaneSurface([curveLoop])

gmsh.model.occ.synchronize()

lxx = [l1,l3,l9,l11]
lyy = [l2,l4,l10,l12]
lzz = [l5,l6,l7,l8]

for i in lxx:
    gmsh.model.mesh.setTransfiniteCurve(i,nx)


for i in lyy:
    gmsh.model.mesh.setTransfiniteCurve(i,ny)
    

for i in lzz:
    gmsh.model.mesh.setTransfiniteCurve(i,nz)

ff = [f1,f2,f3,f4,f5,f6]
for i in ff:
    gmsh.model.mesh.setTransfiniteSurface(i)

gmsh.model.mesh.setTransfiniteVolume(vol)

#gmsh.option.setNumber("Mesh.Algorithm", 11)
gmsh.option().setNumber("Mesh.RecombineAll", 1)
#gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
gmsh.option.setNumber("Mesh.ElementOrder", 1)

gmsh.model.mesh.generate(1)
gmsh.model.mesh.generate(2)
gmsh.model.mesh.generate(3)
#gmsh.model.mesh.refine()
#gmsh.model.mesh.reverse([])


gmsh.fltk.run()

pathname = os.path.dirname(sys.argv[0])
currPath = os.path.abspath(pathname)
fesys= HierAMuS.FEMPy(pathname, "firstvolume")
fesys.setStaticSolutionState()
fesys.setSolver(4)
fesys.getMacroCommands().setLogLevel(fesys.BasicLog(), fesys.BasicLog())


gm = fesys.getMeshCommands().getFromGMESH()
gm.addGeomFromGmsh(gmsh)
fesys.getMeshCommands().getGeometryCommands().checkGeometry()

gm.addBrickVolumeFiniteElements(gmsh,1,vol,1)


fesys.getMeshCommands().getElementFormulations().addEL300_3DSolid(1,meshiddisp=1,disporder=disporder,mode=2)
#fesys.getMeshCommands().getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1,E=E,nu=nu)
fesys.getMeshCommands().getMaterialFormulations().addMA3_SmallStrainPlasticity(1,E=E,nu=nu,y0=fy,yinf=100,xh=xh,xd=0,eta=0)
fesys.getMeshCommands().addMaterial(1,1,1)

fesys.getMeshCommands().setDegreesOfFreedom()

bounfaces = gm.getFaceNumbers(gmsh,f5,4,1)
fesys.getMeshCommands().getBoundaryConditions().singleBC(fesys.getMeshCommands().getGeometryCommands().faceType(),bounfaces,1,[1,1,1],disporder)

loadfaces = gm.getFaceNumbers(gmsh,f3,4,1)
fesys.getMeshCommands().getBoundaryConditions().singleLoad(fesys.getMeshCommands().getGeometryCommands().faceType(),loadfaces,1,[0,1,0],0,shapeorder=disporder)


gmsh.finalize()
fesys.getMacroCommands().sparseSetUp()

fesys.getMacroCommands().setPropFunction(0,lambda t:t,0,0.5)
fesys.getMacroCommands().setPropFunction(0,lambda t:1-t,0.5,1)


fesys.getMacroCommands().setLogLevel(fesys.NoLog(), fesys.NoLog())

steps = 40
fesys.getMacroCommands().setDt(1/steps)
fesys.getPlotCommands().toFile()

fesys.getMacroCommands().assembleSolve()
#fesys.getMacroCommands().computeEigenValues(8,60)

fesys.getMacroCommands().setLogLevel(fesys.FullLog(),fesys.FullLog())
for i in range(steps):
    fesys.getMacroCommands().timeincr()
    fesys.getMacroCommands().newton(refResidual=1e-11)
    fesys.getPlotCommands().toFile()


sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(), p5, 1)

#print(p5, sol)



