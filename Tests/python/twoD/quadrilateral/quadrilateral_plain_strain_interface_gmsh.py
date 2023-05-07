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
    nx=20
    ny=20
    esize = 0.5
    disporder = 4
    reffac = 4
    reffac2 = 8

    # initialize gmsh
    gmsh.initialize()
    gmsh.model.add("t11")


    # defining points of geometry
    p1 = gmsh.model.geo.addPoint(0,0,0,esize)
    p2 = gmsh.model.geo.addPoint(L/2,h1/2,0,esize*reffac)
    p3 = gmsh.model.geo.addPoint(L,h1,0,esize/reffac)
    p4 = gmsh.model.geo.addPoint(L,h1+(h2-h1)/2,0,esize)
    p5 = gmsh.model.geo.addPoint(L,h2,0,esize/reffac)
    p6 = gmsh.model.geo.addPoint(L/2,h1+(h2-h1)/2,0,esize*reffac)
    p7 = gmsh.model.geo.addPoint(0,h1,0,esize/reffac2)
    p8 = gmsh.model.geo.addPoint(0,h1/2,0,esize)

    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p5)
    l5 = gmsh.model.geo.addLine(p5, p6)
    l6 = gmsh.model.geo.addLine(p6, p7)
    l7 = gmsh.model.geo.addLine(p7, p8)
    l8 = gmsh.model.geo.addLine(p8, p1)

    # defining the final face
    cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5, l6, l7, l8])
    pl = gmsh.model.geo.addPlaneSurface([cl])

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.getElements()


    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.option().setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3)

    gmsh.model.mesh.generate(2)
    #gmsh.model.mesh.refine()


    #gmsh.fltk.run()


    # gettting the gmsh mesh handler
    gm = fesys.getMeshCommands().getFromGMESH()

    # adding gmsh mesh as geometry to FEMPy
    gm.addGeomFromGmsh(gmsh)
    fesys.getMeshCommands().getGeometryCommands().checkGeometry()

    fesys.getMacroCommands().setLogLevel(fesys.BasicLog(), fesys.BasicLog())
    fesys.getMacroCommands().printInfo()

    sys.exit(0)

    # adding the finite elements as 2D Quadrilaterls
    eltype = gmsh.model.mesh().getElementType("Quadrangle", 1)
    gm.addFiniteElementsFromTag(gmsh, eltype, pl, 1)

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
    bedges = gm.getEdgeNumbers(gmsh, l7) + gm.getEdgeNumbers(gmsh, l8)
    geo = fesys.getMeshCommands().getGeometryCommands()
    fesys.getMeshCommands().getBoundaryConditions().singleBC(geo.edgeType(), bedges, 1, [1,1,1], disporder)

    bedges = gm.getEdgeNumbers(gmsh, l3) + gm.getEdgeNumbers(gmsh, l4)
    fesys.getMeshCommands().getBoundaryConditions().singleLoad(geo.edgeType(),bedges,1,[0,6.25,0],0,shapeorder=disporder)

    fesys.getMacroCommands().sparseSetUp()

    fesys.getMacroCommands().setPropFunction(0)

    #fesys.getMacroCommands().newton(refResidual=1e-10)

    steps = 20
    fesys.getMacroCommands().setDt(1.0/steps)
    for i in range(steps):
        fesys.getMacroCommands().timeincr()
        fesys.getMacroCommands().newton(maxIteration=50,refResidual=1e-7)

    #fesys.getMacroCommands().printInfo()

    print(p2)
    sol = fesys.getMacroCommands().getSolution(geo.vertexType(), geomNumber=p4, meshId=1)
    print(sol)

    fesys.getPlotCommands().toFile()


firstModel()