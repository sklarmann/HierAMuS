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

nx = 30
ny = 30


p1 = gmsh.model.occ.addPoint(0,0,0)
p2 = gmsh.model.occ.addPoint(L,0,0)
p3 = gmsh.model.occ.addPoint(L,h,0)
p4 = gmsh.model.occ.addPoint(0,h,0)

l1 = gmsh.model.occ.addLine(p1,p2)
l2 = gmsh.model.occ.addLine(p2,p3)
l3 = gmsh.model.occ.addLine(p3,p4)
l4 = gmsh.model.occ.addLine(p4,p1)

cl = gmsh.model.occ.addCurveLoop([l1,l2,l3,l4])

f1 = gmsh.model.occ.addPlaneSurface([cl])

gmsh.model.occ.synchronize()

gmsh.option.setNumber("Mesh.Algorithm",8)
gmsh.option.setNumber("Mesh.RecombineAll", 1)
#gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)

gmsh.model.mesh.setTransfiniteCurve(l1,nx)
gmsh.model.mesh.setTransfiniteCurve(l2,ny)
gmsh.model.mesh.setTransfiniteCurve(l3,nx)
gmsh.model.mesh.setTransfiniteCurve(l4,ny)
gmsh.model.mesh.setTransfiniteSurface(f1)

gmsh.model.mesh.generate(2)


pathname = os.path.dirname(sys.argv[0])
currPath = os.path.abspath(pathname)
fesys= HierAMuS.FEMPy(pathname, "firstvolume")
fesys.setSolutionState()
fesys.setSolver(3)
fesys.getMacroCommands().setLogLevel(fesys.FullLog(), fesys.FullLog())

gm = fesys.getMeshCommands().getFromGMESH()
geo = fesys.getMeshCommands().getGeometryCommands()


gm.addGeomFromGmsh(gmsh)
geo.checkGeometry()

gm.addQuadrilateralFiniteElements(gmsh,1,f1,1)

fesys.getMeshCommands().getElementFormulations().addEL201_2DShell(1,1,1,1)
fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(1,100,0.3,1,1)
fesys.getMeshCommands().addMaterial(1,1,1)

fesys.getMeshCommands().setDegreesOfFreedom()


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

