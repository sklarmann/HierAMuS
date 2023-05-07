# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from operator import ge
import sys, os
import HierAMuS
import numpy as np

E=100
nu=0.3


disporder = 1

pathname = os.path.dirname(sys.argv[0])
currPath = os.path.abspath(pathname)
fesys= HierAMuS.FEMPy(pathname, "firstvolume")
fesys.setSolutionState()
fesys.setSolver(3)
fesys.getMacroCommands().setLogLevel(fesys.FullLog(), fesys.FullLog())


geo = fesys.getMeshCommands().getGeometryCommands()

geo.addVertex(1,0,0,0)
geo.addVertex(2,1,0,0)
geo.addVertex(3,1,1,0)
geo.addVertex(4,0,1,0)
geo.addVertex(5,2,0,0)
geo.addVertex(6,2,1,0)

geo.addLinearEdgeGeo(5,[1,4])
geo.addLinearEdgeGeo(6,[5,6])

geo.addQuadrilateralFace(1,[1,2,3,4])
geo.addQuadrilateralFace(2,[2,5,6,3])

fesys.getMeshCommands().getGeometryCommands().checkGeometry()

fem = fesys.getMeshCommands().getElementCommands()
fem.addFace(1,[1,2])



#fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(number=1, E=E, nu=nu, thickness=1, plainstrain=0)

fesys.getMeshCommands().getMaterialFormulations().addMA1_2D_PlainStrain_3D(1, 2)
Lambda = E*nu/((1+nu)*(1-2*nu))
fesys.getMeshCommands().getMaterialFormulations().addMA2_3D_NeoHook(2, Lambda, E/2/(1+nu))
#fesys.getMeshCommands().getMaterialFormulations().addMA3_SmallStrainPlasticity(2, E=E, nu=nu, y0=E/100, yinf=0, xh=E/10, xd=0, eta=0)

fesys.getMeshCommands().getElementFormulations().addEL201_2DShell(num=1, meshiddisp=1, disporder=disporder, mode=3)
fesys.getMeshCommands().addMaterial(matNum=1, matFormNum=1, elemFormNum=1)

fesys.getMeshCommands().setDegreesOfFreedom()


fesys.getMeshCommands().getBoundaryConditions().singleBC(fesys.getMeshCommands().getGeometryCommands().edgeType(), number=5, meshId=1, dofs=[1,1,1], shapeOrder=disporder)
fesys.getMacroCommands().printInfo()
fesys.getMeshCommands().getBoundaryConditions().singleLoad(fesys.getMeshCommands().getGeometryCommands().edgeType(), number=6, meshId=1, load=[0,20,0], propnum=0,shapeorder=disporder)

fesys.getMacroCommands().sparseSetUp()
fesys.getMacroCommands().setPropFunction(0)


#fesys.getMacroCommands().setLogLevel(fesys.Basi(), fesys.NoLog())

steps = 10
fesys.getMacroCommands().setDt(1/steps)
fesys.getPlotCommands().toFile()

fesys.getMacroCommands().assembleSolve()
#fesys.getMacroCommands().computeEigenValues(8,60)
for i in range(steps):
    fesys.getMacroCommands().timeincr()
    fesys.getMacroCommands().newton(refResidual=1e-11)
    fesys.getPlotCommands().toFile()


#sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(), p5, 1)





