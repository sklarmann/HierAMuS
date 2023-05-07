# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import gmsh
import sys, os
import HierAMuS
import numpy as np

gmsh.initialize()

gmsh.model.add("t11")


b=5
h=5
r=1

disporder=1

E=100
nu=0.49
e=1

p1 = gmsh.model.occ.addPoint(0, 0, 0,e)
p2 = gmsh.model.occ.addPoint(1, 0, 0,e/20)
p3 = gmsh.model.occ.addPoint(5, 0, 0,e)
p4 = gmsh.model.occ.addPoint(5, 5, 0,e)
p5 = gmsh.model.occ.addPoint(0, 5, 0,e)
p6 = gmsh.model.occ.addPoint(0, 1, 0,e/10)

l1=gmsh.model.occ.addLine(p6, p5)
l2=gmsh.model.occ.addLine(p5, p4)
l3=gmsh.model.occ.addLine(p3, p4)
l4=gmsh.model.occ.addLine(p3, p2)
circ = gmsh.model.occ.addCircleArc(p2, p1, p6)

curveLoop = gmsh.model.occ.addCurveLoop([-l1, -l2, l3, -l4, -circ])
surface = gmsh.model.occ.addPlaneSurface([curveLoop])

gmsh.model.occ.synchronize()


gmsh.option.setNumber("Mesh.Algorithm", 11)
gmsh.option().setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
gmsh.option.setNumber("Mesh.ElementOrder", 1)

gmsh.model.mesh.generate(2)

gmsh.model.mesh.reverse([])

gmsh.fltk().run()



pathname = os.path.dirname(sys.argv[0])
currPath = os.path.abspath(pathname)
fesys= HierAMuS.FEMPy(pathname, "firstvolume")
fesys.setStaticHomogenizationSolutionState()
fesys.setSolver(3)
fesys.getMacroCommands().setLogLevel(fesys.BasicLog(), fesys.BasicLog())


gm = fesys.getMeshCommands().getFromGMESH()
gm.addGeomFromGmsh(gmsh)
fesys.getMeshCommands().getGeometryCommands().checkGeometry()

gm.addQuadrilateralFiniteElements(gmsh, 1, p1, 1)

#fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(number=1, E=E, nu=nu, thickness=1, plainstrain=0)

fesys.getMeshCommands().getMaterialFormulations().addMA1_2D_PlainStrain_3D(1, 2)
Lambda = E*nu/((1+nu)*(1-2*nu))
fesys.getMeshCommands().getMaterialFormulations().addMA2_3D_NeoHook(2, Lambda, E/2/(1+nu))
#fesys.getMeshCommands().getMaterialFormulations().addMA3_SmallStrainPlasticity(2, E=E, nu=nu, y0=E/100, yinf=0, xh=E/10, xd=0, eta=0)

fesys.getMeshCommands().getElementFormulations().addEL201_2DShell(num=1, meshiddisp=1, disporder=disporder, mode=3)
fesys.getMeshCommands().addMaterial(matNum=1, matFormNum=1, elemFormNum=1)

fesys.getMeshCommands().setDegreesOfFreedom()

enums = gm.getEdgeNumbers(gmsh, tag=l1, order=1)
fesys.getMeshCommands().getBoundaryConditions().singleBC(fesys.getMeshCommands().getGeometryCommands().edgeType(), number=enums, meshId=1, dofs=[1,0,1], shapeOrder=disporder)
enums = gm.getEdgeNumbers(gmsh, tag=l4, order=1)
fesys.getMeshCommands().getBoundaryConditions().singleBC(fesys.getMeshCommands().getGeometryCommands().edgeType(), number=enums, meshId=1, dofs=[0,1,1], shapeOrder=disporder)

enums = gm.getEdgeNumbers(gmsh, tag=l2, order=1)
fesys.getMeshCommands().getBoundaryConditions().singleLoad(fesys.getMeshCommands().getGeometryCommands().edgeType(), number=enums, meshId=1, load=[0,20,0], propnum=0,shapeorder=disporder)
gmsh.finalize()
fesys.getMacroCommands().sparseSetUp()
fesys.getMacroCommands().setPropFunction(0)


fesys.getMacroCommands().setLogLevel(fesys.NoLog(), fesys.NoLog())

steps = 5
fesys.getMacroCommands().setDt(1/steps)
fesys.getPlotCommands().toFile()

fesys.getMacroCommands().assembleSolve()
fesys.getMacroCommands().computeEigenValues(8,60)
for i in range(steps):
    fesys.getMacroCommands().timeincr()
    fesys.getMacroCommands().newton(refResidual=1e-11)
    fesys.getPlotCommands().toFile()


sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(), p5, 1)

print(p5, sol)



