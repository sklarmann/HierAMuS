# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from operator import ge
import sys, os
import HierAMuS
import numpy as np
import gmsh


gmsh.initialize()
gmsh.model.add('test')

L=1
h=1
b=1

nx = 3
ny = 3
nz = 3


p1 = gmsh.model.occ.addPoint(0,0,0)
p2 = gmsh.model.occ.addPoint(L,0,0)
p3 = gmsh.model.occ.addPoint(L,h,0)
p4 = gmsh.model.occ.addPoint(0,h,0)
p5 = gmsh.model.occ.addPoint(0,0,b)
p6 = gmsh.model.occ.addPoint(L,0,b)
p7 = gmsh.model.occ.addPoint(L,h,b)
p8 = gmsh.model.occ.addPoint(0,h,b)

pl = [p1,p2,p3,p4,p5,p6,p7,p8]

l1 = gmsh.model.occ.addLine(p1,p2)
l2 = gmsh.model.occ.addLine(p2,p3)
l3 = gmsh.model.occ.addLine(p3,p4)
l4 = gmsh.model.occ.addLine(p4,p1)

l5 = gmsh.model.occ.addLine(p1,p5)
l6 = gmsh.model.occ.addLine(p2,p6)
l7 = gmsh.model.occ.addLine(p3,p7)
l8 = gmsh.model.occ.addLine(p4,p8)


l9  = gmsh.model.occ.addLine(p5,p6)
l10 = gmsh.model.occ.addLine(p6,p7)
l11 = gmsh.model.occ.addLine(p7,p8)
l12 = gmsh.model.occ.addLine(p8,p5)


llx = [l1,l3,l9,l11]
lly = [l2,l4,l10,l12]
llz = [l5,l6,l7,l8]

cl1 = gmsh.model.occ.addCurveLoop([l1,l2,l3,l4])
cl2 = gmsh.model.occ.addCurveLoop([l1,l6,-l9,-l5])
cl3 = gmsh.model.occ.addCurveLoop([l2,l7,-l10,-l6])
cl4 = gmsh.model.occ.addCurveLoop([l3,l8,-l11,-l7])
cl5 = gmsh.model.occ.addCurveLoop([l4,l5,-l12,-l8])
cl6 = gmsh.model.occ.addCurveLoop([l9,l10,-l11,-l12])

f1 = gmsh.model.occ.addPlaneSurface([cl1])
f2 = gmsh.model.occ.addPlaneSurface([cl2])
f3 = gmsh.model.occ.addPlaneSurface([cl3])
f4 = gmsh.model.occ.addPlaneSurface([cl4])
f5 = gmsh.model.occ.addPlaneSurface([cl5])
f6 = gmsh.model.occ.addPlaneSurface([cl6])

faces = [f1,f2,f3,f4,f5,f6]

vl = gmsh.model.occ.addSurfaceLoop([f1,f2,f3,f4,f5,f6])
gmsh.model.occ.addVolume([vl])

gmsh.model.occ.synchronize()


gmsh.option.setNumber("Mesh.Algorithm",8)
gmsh.option.setNumber("Mesh.RecombineAll", 1)

for i in llx:
    gmsh.model.mesh.setTransfiniteCurve(i,nx)
for i in lly:
    gmsh.model.mesh.setTransfiniteCurve(i,ny)
for i in llz:
    gmsh.model.mesh.setTransfiniteCurve(i,nz)
for i in faces:
    gmsh.model.mesh.setTransfiniteSurface(i)
    
gmsh.model.mesh.setTransfiniteVolume(vl)

gmsh.model.mesh.generate(3)






pathname = os.path.dirname(sys.argv[0])
currPath = os.path.abspath(pathname)
fesys= HierAMuS.FEMPy(pathname, "firstvolume")
fesys.setSolutionState()
fesys.setSolver(3)
fesys.getMacroCommands().setLogLevel(fesys.BasicLog(), fesys.BasicLog())

gm = fesys.getMeshCommands().getFromGMESH()
geo = fesys.getMeshCommands().getGeometryCommands()



gm.addGeomFromGmsh(gmsh)
geo.checkGeometry()

#fesys.getMacroCommands().printInfo()

#gmsh.fltk.run()

gm.addBrickVolumeFiniteElements(gmsh, 1, f1, 1)

fesys.getMeshCommands().getElementFormulations().addEL300_3DSolid(1,1,3,1)
fesys.getMeshCommands().getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1,100,0.3)
fesys.getMeshCommands().addMaterial(1,1,1)





fesys.getMeshCommands().setDegreesOfFreedom()

#fesys.getMacroCommands().printInfo()


fnums = gm.getFaceNumbers(gmsh,f5,4,1)
fesys.getMeshCommands().getBoundaryConditions().singleBC(geo.faceType(), fnums, 1, [1,1,1], 1)

fnums = gm.getFaceNumbers(gmsh,f3,4,1)
fesys.getMeshCommands().getBoundaryConditions().singleLoad(geo.faceType(), fnums, 1, [0,1,0], 1, shapeorder=1)

fesys.getMacroCommands().sparseSetUp()

fesys.getMacroCommands().setPropFunction(1)
fesys.getMacroCommands().setDt(1)
fesys.getMacroCommands().timeincr()

fesys.getMacroCommands().assembleSolve()

fesys.getPlotCommands().toFile()
#fesys.getMacroCommands().printInfo()

sys.exit(0)



et = gm.getEdgeNumbers(gmsh,l4,1)
print(et)

fesys.getMeshCommands().getBoundaryConditions().singleBC(geo.edgeType(),et,1,[1,1,1],1)
et = gm.getEdgeNumbers(gmsh,l2,1)
fesys.getMeshCommands().getBoundaryConditions().singleLoad(geo.edgeType(),et,1,[0,1,0],1)

fesys.getMacroCommands().sparseSetUp()

fesys.getMacroCommands().setPropFunction(1)
fesys.getMacroCommands().setDt(1)
fesys.getMacroCommands().timeincr()

fesys.getMacroCommands().newton()
#fesys.getMacroCommands().printInfo()

