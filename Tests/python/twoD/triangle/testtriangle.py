# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

# from operator import ge
import os
import sys

import HierAMuS
import numpy as np

Emod=1
nue=0.0


disporder = 2

pathname = os.path.dirname(sys.argv[0])
currPath = os.path.abspath(pathname)
fesys= HierAMuS.FEMPy(pathname, "firstvolume1")
fesys.setStaticSolutionState()
fesys.getMacroCommands().setLogLevel(fesys.BasicLog(), fesys.FullLog())


geo = fesys.getMeshCommands().getGeometryCommands()

geo.addVertex(1,0,0,0)
geo.addVertex(2,1,0,0)
geo.addVertex(3,0,1,0)
geo.addVertex(4,1,1,0)

# geo.addVertex(5,1,2,0)
# geo.addVertex(6,0,2,0)


#geo.addLinearEdgeGeo(1,[1,2])
#geo.addLinearEdgeGeo(2,[2,3])
#geo.addLinearEdgeGeo(3,[3,1])
# geo.addLinearEdgeGeo(4,[2,4])
# geo.addLinearEdgeGeo(5,[4,3])

geo.addTriangleFace(1,[1,2,3])
geo.addTriangleFace(2,[4,3,2])
# geo.addQuadrilateralFace(1,[1,2,4,3])
# geo.addQuadrilateralFace(2,[4,3,5,6])

# geo.addQuadrilateralFace(1,[1,2,4,3])
# geo.addQuadrilateralFace(2,[3,5,6,4])

# geo.addQuadrilateralFace(1,[1,2,4,3])
# geo.addQuadrilateralFace(2,[5,6,4,3])

# geo.addQuadrilateralFace(1,[1,2,4,3])
# geo.addQuadrilateralFace(2,[6,4,3,5])


fesys.getMeshCommands().getGeometryCommands().checkGeometry()

fem = fesys.getMeshCommands().getElementCommands()
fem.addFace(1,[1,2])


#fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(number=1, E=E, nu=nu, thickness=1, plainstrain=0)

fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(1,E=Emod,nu=nue,thickness=1,plainstrain=0)
# Lambda = E*nu/((1+nu)*(1-2*nu))
# fesys.getMeshCommands().getMaterialFormulations().addMA2_3D_NeoHook(2, Lambda, E/2/(1+nu))
#fesys.getMeshCommands().getMaterialFormulations().addMA3_SmallStrainPlasticity(2, E=E, nu=nu, y0=E/100, yinf=0, xh=E/10, xd=0, eta=0)

# fesys.getMeshCommands().getElementFormulations().addEL201_2DShell(num=1, meshiddisp=1, disporder=disporder, mode=1)
fesys.getMeshCommands().getElementFormulations().addEL205_HDivTest(1,100,disporder=disporder,stressorder=disporder-1,mode=2,meshiddisp=1,meshidstress=2,E=1,nu=0.3)
fesys.getMeshCommands().addMaterial(matNum=1, matFormNum=1, elemFormNum=1)

fesys.getMeshCommands().setDegreesOfFreedom()

fesys.getMeshCommands().getBoundaryConditions().singleBC(fesys.getMeshCommands().getGeometryCommands().vertexType(), number=[1], meshId=1, dofs=[1,1,1], shapeOrder=disporder)
fesys.getMeshCommands().getBoundaryConditions().singleBC(fesys.getMeshCommands().getGeometryCommands().vertexType(), number=[2], meshId=1, dofs=[0,1,1], shapeOrder=disporder)
fesys.getMeshCommands().getBoundaryConditions().singleLoad(fesys.getMeshCommands().getGeometryCommands().vertexType(), number=[3,4], meshId=1, load=[0,1,0], propnum=0,shapeorder=disporder)
# fesys.getMacroCommands().printInfo()

fesys.getMacroCommands().sparseSetUp()
fesys.getMacroCommands().setPropFunction(0)


#fesys.getMacroCommands().setLogLevel(fesys.Basi(), fesys.NoLog())

steps = 1
fesys.getMacroCommands().setDt(1/steps)
fesys.getMacroCommands().timeincr()
fesys.getMacroCommands().assembleSolve()

# fesys.getPlotCommands().toFile()
# fesys.getMacroCommands().printSpMat()
sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(),geomNumber=3,meshId=1)
print(sol)
sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().edgeType(),geomNumber=2,meshId=1)
print(sol)
#fesys.getMacroCommands().computeEigenValues(8,60)
# for i in range(steps):
#     fesys.getMacroCommands().timeincr()
    # fesys.getMacroCommands().newton(refResidual=1e-11)
    # fesys.getPlotCommands().toFile()


#sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(), p5, 1)





