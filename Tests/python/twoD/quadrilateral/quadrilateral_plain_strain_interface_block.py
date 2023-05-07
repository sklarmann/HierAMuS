# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import gmsh
import sys, os
import HierAMuS


def firstModel():
    pathname = os.path.dirname(sys.argv[0])
    currPath = os.path.abspath(pathname)


    fesys = HierAMuS.FEMPy(currPath, 'firsttest')
    fesys.setStaticSolutionState()
    fesys.setSolver(2)
    fesys.getMacroCommands().setLogLevel(fesys.FullLog(), fesys.FullLog())

    # geometry and element parameters

    L=48
    h1=44
    h2=60
    nx=200
    ny=200
    disporder = 2

    geoB = fesys.getGeomBuilder()

    # Lower Block
    # Vertices definieren
    v1 = geoB.getVertex()
    v1.setCoordinates(0,0,0)
    v2 = geoB.getVertex()
    v2.setCoordinates(L,h1,0)
    v3 = geoB.getVertex()
    v3.setCoordinates(L,h2,0)
    v4 = geoB.getVertex()
    v4.setCoordinates(0,h1,0)
    # Edges festlegen aus Vertices.
    e1 = geoB.getEdge()
    e1.setStartEnd(v2,v1)
    e2 = geoB.getEdge()
    e2.setStartEnd(v3,v2)
    e3 = geoB.getEdge()
    e3.setStartEnd(v4,v3)
    e4 = geoB.getEdge()
    e4.setStartEnd(v4,v1)
    q1 = geoB.getQuad()
    q1.setVertsEdges([v1,v2,v3,v4],[e1,e2,e3,e4])
    q1.setDivision(nx,ny)
    geoB.process()

    fesys.getMeshCommands().getElementCommands().addFace(1,q1.getQuadList())
    meshCMD = fesys.getMeshCommands()

    # setting up the element and material formulation
    fesys.getMeshCommands().getElementFormulations().addEL201_2DShell(num=1, meshiddisp=1, disporder=disporder, mode=3)
    fesys.getMeshCommands().getMaterialFormulations().addMA1_2D_PlainStrain_3D(1, 2)
    #fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(number=1, E=10000, nu=0.3, thickness=1, plainstrain=1)
    #fesys.getMeshCommands().getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(number=2, E=10000, nu=0.3)
    E=240.5654
    v=0.4999
    G = E/2/(1+v)
    lam = E/3/(1-2*v)
    fesys.getMeshCommands().getMaterialFormulations().addMA2_3D_NeoHook(number=2, Lambda=lam, G=G)
    fesys.getMeshCommands().addMaterial(matNum=1, matFormNum=1, elemFormNum=1)

    fesys.getMeshCommands().setDegreesOfFreedom()

    meshCMD.getBoundaryConditions().singleBC(meshCMD.getGeometryCommands().edgeType(),number=e4.EdgeList,meshId=1,dofs=[1,1,1],shapeOrder=disporder)
    meshCMD.getBoundaryConditions().singleLoad(meshCMD.getGeometryCommands().edgeType(),number=e2.EdgeList,meshId=1,load=[0,6.25,0],shapeorder=disporder,propnum=0)





    fesys.getMacroCommands().sparseSetUp()

    fesys.getMacroCommands().setPropFunction(0)

    #fesys.getMacroCommands().newton(refResidual=1e-10)

    steps = 20
    fesys.getMacroCommands().setDt(1.0/steps)
    for i in range(steps):
        fesys.getMacroCommands().timeincr()
        fesys.getMacroCommands().newton(maxIteration=50,refResidual=1e-7)

    #fesys.getMacroCommands().printInfo()
    macroCMD = fesys.getMacroCommands()

    sol = macroCMD.getSolution(meshCMD.getGeometryCommands().vertexType(),geomNumber=v3.num,meshId=1)
    print(sol)

    fesys.getPlotCommands().toFile()


firstModel()