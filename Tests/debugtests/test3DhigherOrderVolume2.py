# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import HierAMuS

import os
import sys
import time
import gmsh
import math
import matplotlib.pyplot as plt
# time.sleep(20)

x = []
y = []

for i in range(1):
    L=10
    b=1
    h=1

    nx = 1
    ny = 1
    order = 1

    # Ermittle Pfad des Skripts
    print('sys.argv[0] =', sys.argv[0])
    pathname = os.path.dirname(sys.argv[0])
    print('path =', pathname)
    print('full path =', os.path.abspath(pathname))

    currPath = os.path.abspath(pathname)

    fesys = HierAMuS.FEMPy(currPath, "test3DBeamInterface")
    fesys.setStaticSolutionState()
    fesys.setSolver(1)
    Macro = fesys.getMacroCommands()
    Macro.setLogLevel(fesys.FullLog(), fesys.NoLog())


    Mesh = fesys.getMeshCommands()
    Geo = Mesh.getGeometryCommands()

    Geo.addVertex(1,0,0,0)
    Geo.addVertex(2,L,0,0)
    Geo.addVertex(3,L,b,0)
    Geo.addVertex(4,0,b,0)
    Geo.addVertex(5,0,0,h)
    Geo.addVertex(6,L,0,h)
    Geo.addVertex(7,L,b,h)
    Geo.addVertex(8,0,b,h)

    Geo.addLinearEdgeGeo(1,[1,2])
    Geo.addLinearEdgeGeo(2,[2,3])
    Geo.addLinearEdgeGeo(3,[3,4])
    Geo.addLinearEdgeGeo(4,[1,4])


    Geo.addLinearEdgeGeo(5,[1,5])
    Geo.addLinearEdgeGeo(6,[2,6])
    Geo.addLinearEdgeGeo(7,[3,7])
    Geo.addLinearEdgeGeo(8,[4,8])

    Geo.addLinearEdgeGeo(9, [6,5])
    Geo.addLinearEdgeGeo(10,[6,7])
    Geo.addLinearEdgeGeo(11,[7,8])
    Geo.addLinearEdgeGeo(12,[8,5])

    Geo.addQuadrilateralFace(1,[1,2,3,4],[1,2,3,4])
    Geo.addQuadrilateralFace(2,[1,2,6,5],[1,6,9,5])
    Geo.addQuadrilateralFace(3,[2,3,7,6],[2,7,10,6])
    Geo.addQuadrilateralFace(4,[3,4,8,7],[3,7,11,7])
    Geo.addQuadrilateralFace(5,[1,5,8,4],[5,12,8,4])
    Geo.addQuadrilateralFace(6,[6,7,8,5],[10,11,12,9])

    Geo.addLinearBrick(1,[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8,9,10,11,12],[1,2,3,4,5,6])

    Mesh.getElementCommands().addVolume(1,1)



    #Macro.printInfo()


    # Adding the Volume Elements and assign the material
    Mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1, E=1000, nu=0.0)
    Mesh.getElementFormulations().addEL300_3DSolid(1, meshiddisp=1, disporder=order, mode=1)
    Mesh.addMaterial(1, 1, 1)

    # Adding the interface element on the load side and assign the material
    Mesh.setDegreesOfFreedom()

    Mesh.getBoundaryConditions().singleBC(Geo.faceType(), number=5, meshId=1, dofs=[1,1,1], shapeOrder=order)
    Mesh.getBoundaryConditions().singleLoad(Geo.faceType(), number=3, meshId=1, load=[0,1,0], propnum=1, shapeorder=order)
    #Mesh.getBoundaryConditions().singleLoad(Geo.vertexType(), number=vnum, meshId=2, load=[1,0,0], propnum=1)

    Macro.sparseSetUp()

    Macro.setPropFunction(1)
    Macro.setPropFunction(0)
    Macro.setDt(1)
    Macro.timeincr()

    Macro.setLogLevel(fesys.FullLog(),fesys.FullLog())

    #Macro.newton(maxIteration=2)
    #Macro.printInfo()
    Macro.assembleSolve()
    Macro.printInfo()

    fesys.getPlotCommands().toFile()

    sol=Macro.getSolution(Geo.vertexType(),2,1)
    x.append(i+1)
    y.append(sol[1])
    #print(sol)
    #sol=Macro.getSolution(Geo.vertexType(),vnum,2)
    #print(sol)

print(y)
fig, ax = plt.subplots()
ax.plot(x, y, linewidth=2.0)
#plt.show()