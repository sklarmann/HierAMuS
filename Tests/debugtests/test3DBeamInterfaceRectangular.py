# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import gmsh

import os, sys

import HierAMuS


def run():
    gmsh.initialize()

    nlT = 6
    ncT = 8
    nb = 6
    nl = 40

    ll = 0.01
    L = 10
    L-=2*ll

    disporder = 2
    warptype = 3

    gmsh.model.add("RectangularBeam")

    gmsh.model.occ.addPoint(ll, -0.5, -0.5)
    gmsh.model.occ.addPoint(ll,  0.5, -0.5)
    gmsh.model.occ.addPoint(ll,  0.5, 0.5)
    gmsh.model.occ.addPoint(ll, -0.5, 0.5)

    ps = gmsh.model.occ.addPoint(0, 0, 0)
    pe = gmsh.model.occ.addPoint(L+2*ll, 0, 0)

    gmsh.model.occ.addLine(1, 2)
    gmsh.model.occ.addLine(2, 3)
    gmsh.model.occ.addLine(3, 4)
    gmsh.model.occ.addLine(4, 1)



    lines = [1,2,3,4]
    gmsh.model.occ.addCurveLoop(lines)

    gmsh.model.occ.addPlaneSurface([1])
    

    backfaces = [1]


    gmsh.model.occ.extrude([(2,1)], L, 0, 0,[nl],recombine=True)
    gmsh.model.occ.synchronize()


    for i in lines:
        gmsh.model.mesh.setTransfiniteCurve(i, nlT)
    for i in backfaces:
        gmsh.model.mesh.setTransfiniteSurface(i)


    frontfaces = [6]
    volumes = [1]

    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)

    gmsh.model.mesh.generate(0)
    gmsh.model.mesh.generate(3)

    #gmsh.model.mesh.reverse([(2,1),(2,3)])
    gmsh.model.mesh.reverse([(2,2)])
    gmsh.model.mesh.removeDuplicateNodes()


    nt, nc, parm = gmsh.model.mesh.getNodes(dim=0,tag = pe)
    print(nt)

    gmsh.model.occ.synchronize()

    #gmsh.fltk().run()


    pathname = os.path.dirname(sys.argv[0])
    currPath = os.path.abspath(pathname)
    fesys= HierAMuS.FEMPy(pathname, "RectangularCantilever")
    fesys.setStaticSolutionState()
    fesys.setSolver(4)

    mesh = fesys.getMeshCommands()
    gm = mesh.getFromGMESH()
    geo = mesh.getGeometryCommands()
    macro = fesys.getMacroCommands()

    macro.setLogLevel(fesys.FullLog(),fesys.FullLog())

    gm.addGeomFromGmsh(gmsh)
    geo.checkGeometry()

    gm.addBrickVolumeFiniteElements(gmsh, 1, [1], 1)

    mesh.getElementFormulations().addEL300_3DSolid(1,1,disporder,1)
    mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1,100,0.0)
    mesh.addMaterial(1,1,1)



    # Back Face
    nt,cn,t =gmsh.model.mesh.getNodes(dim=0,tag=ps)
    vnum1 = nt[0]
    fnums = gm.getFaceNumbers(gmsh,backfaces[0],4,1)
    backmat=[]
    for i in fnums:
        backmat.append(1)
    mesh.getElementCommands().add3DInterfaceElement(materialNumber=2, specialNumber=1, vertNum=vnum1, faceList=fnums,faceMaterialList=backmat)
    mesh.getElementFormulations().addEL302_BeamCoupling3D(num=2, disporder=disporder, meshiddisp=1, meshidrot=2,mode=3,warptype=warptype)
    mesh.addMaterial(2,1,2)

    # Front Face
    nt,cn,t =gmsh.model.mesh.getNodes(dim=0,tag=pe)
    vnum2 = nt[0]
    fnums = gm.getFaceNumbers(gmsh,frontfaces[0],4,1)
    backmat = []
    for i in fnums:
        backmat.append(1)
    mesh.getElementCommands().add3DInterfaceElement(materialNumber=3, specialNumber=2, vertNum=vnum2, faceList=fnums,faceMaterialList=backmat)
    mesh.getElementFormulations().addEL302_BeamCoupling3D(num=3, disporder=disporder, meshiddisp=1, meshidrot=2,mode=3,warptype=warptype)
    mesh.addMaterial(3,1,3)


    mesh.setDegreesOfFreedom()

    mesh.getBoundaryConditions().singleBC(geo.vertexType(), vnum1, 1, [1,1,1], shapeOrder=1)
    mesh.getBoundaryConditions().singleBC(geo.vertexType(), vnum1, 2, [1,1,1], shapeOrder=1)

    mesh.getBoundaryConditions().singleLoad(geo.vertexType(), vnum2, 1, [0,1,0], 1)

    macro.sparseSetUp()

    macro.setPropFunction(1)
    macro.setDt(1)
    macro.timeincr()

    macro.assembleSolve()

    fesys.getPlotCommands().toFile()



run()