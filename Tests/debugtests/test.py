# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import os,sys
import HierAMuS


pathname = os.path.dirname(sys.argv[0])
currPath = os.path.abspath(pathname)

meshid = 1
order =1
L  = 10
h1 = 1
nx = 2
ny = 3
fesys = HierAMuS.FEMPy(currPath, 'firsttest')
fesys.getMacroCommands().setLogLevel(fesys.FullLog(),fesys.FullLog())
fesys.setStaticSolutionState()
fesys.setSolver(1)

geo = fesys.getMeshCommands().getGeometryCommands()

geo.addVertex(1,0,0,0)
geo.addVertex(2,1,0,0)
geo.addVertex(3,1,1,0)
geo.addVertex(4,0,1,0)
geo.addVertex(5,2,0,0)
geo.addVertex(6,2,1,0)

geo.addLinearEdgeGeo(1,[1,4])
geo.addLinearEdgeGeo(2,[5,6])
geo.addLinearEdgeGeo(3,[2,3])

geo.addQuadrilateralFace(1,[1,2,3,4])
geo.addQuadrilateralFace(2,[2,5,6,3])

geo.checkGeometry()

elem = fesys.getMeshCommands().getElementCommands()
elem.addFace(1,[1,2])

fesys.getMeshCommands().getElementFormulations().addEL201_2DShell(num=1,meshiddisp=1,disporder=2,mode=1)
#fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(number=1,E=100,nu=0.3,thickness=1,plainstrain=1)
fesys.getMeshCommands().getMaterialFormulations().addMA1_2D_PlainStrain_3D(1,2)
fesys.getMeshCommands().getMaterialFormulations().addMA3_SmallStrainPlasticity(2,100,0.3,1,100,50,0,0)
fesys.getMeshCommands().addMaterial(1,1,1)

fesys.getMeshCommands().setDegreesOfFreedom()

boun = fesys.getMeshCommands().getBoundaryConditions()

boun.singleBC(geo.edgeType(),1,1,[1,1,1],1)
boun.singleLoad(geo.edgeType(),6,1,[0,1,0],1)

fesys.getMacroCommands().sparseSetUp()

fesys.getMacroCommands().setPropFunction(1,lambda t:t,0,5)
fesys.getMacroCommands().setPropFunction(1,lambda t:10-t,5,10)
fesys.getMacroCommands().setDt(1)

for i in range(10):
    fesys.getMacroCommands().timeincr()
    fesys.getMacroCommands().newton()
    fesys.getPlotCommands().toFile()



#fesys.getMacroCommands().printInfo()