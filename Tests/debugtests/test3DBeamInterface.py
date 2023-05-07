# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import gmsh

import os, sys

import HierAMuS

def run():
    gmsh.initialize()

    # Geometry
    tl = 1
    sl = 1
    t = 0.1
    s = 0.1
    ll=0.1
    L=4-ll

    #Mesh
    nthick = 3
    nflansch = 5
    nsteg = 4
    nl = 80

    disporder = 2
    warptype = 3

    x0 = ll



    y1 = -tl/2
    y2 = y1 + tl/2-s/2
    y3 = y2 + s
    y4 = y3 + tl/2 - s/2

    z1 = sl/2 + t
    z2 = z1 - t
    z3 = z2 - sl
    z4 = z3 - t

    gmsh.model.add('IProfBeam')
    gmsh.model.geo.addPoint(x0,y1,z1)
    gmsh.model.geo.addPoint(x0,y2,z1)
    gmsh.model.geo.addPoint(x0,y3,z1)
    gmsh.model.geo.addPoint(x0,y4,z1)

    gmsh.model.geo.addPoint(x0,y1,z2)
    gmsh.model.geo.addPoint(x0,y2,z2)
    gmsh.model.geo.addPoint(x0,y3,z2)
    gmsh.model.geo.addPoint(x0,y4,z2)

    gmsh.model.geo.addPoint(x0,y1,z3)
    gmsh.model.geo.addPoint(x0,y2,z3)
    gmsh.model.geo.addPoint(x0,y3,z3)
    gmsh.model.geo.addPoint(x0,y4,z3)

    gmsh.model.geo.addPoint(x0,y1,z4)
    gmsh.model.geo.addPoint(x0,y2,z4)
    gmsh.model.geo.addPoint(x0,y3,z4)
    gmsh.model.geo.addPoint(x0,y4,z4)

    p1 = gmsh.model.geo.addPoint(L+2*ll,0,0)
    p2 = gmsh.model.geo.addPoint(0,0,0)

    gmsh.model.geo.addLine(1,2)
    gmsh.model.geo.addLine(2,3)
    gmsh.model.geo.addLine(3,4)
    gmsh.model.geo.addLine(5,1)
    gmsh.model.geo.addLine(6,5)
    gmsh.model.geo.addLine(6,2)
    gmsh.model.geo.addLine(6,7)
    gmsh.model.geo.addLine(8,7)
    gmsh.model.geo.addLine(3,7)
    gmsh.model.geo.addLine(8,4)
    gmsh.model.geo.addLine(6,10)
    gmsh.model.geo.addLine(7,11)
    gmsh.model.geo.addLine(9,10)
    gmsh.model.geo.addLine(9,13)
    gmsh.model.geo.addLine(14,10)
    gmsh.model.geo.addLine(11,10)
    gmsh.model.geo.addLine(13,14)
    gmsh.model.geo.addLine(15,14)
    gmsh.model.geo.addLine(15,11)
    gmsh.model.geo.addLine(12,11)
    gmsh.model.geo.addLine(12,16)
    gmsh.model.geo.addLine(15,16)

    cl = gmsh.model.geo.addCurveLoop([-21,20,-19,22])
    gmsh.model.geo.addPlaneSurface([cl])

    cl = gmsh.model.geo.addCurveLoop([-5,6,-1,-4])
    gmsh.model.geo.addPlaneSurface([cl])

    cl = gmsh.model.geo.addCurveLoop([7,-9,-2,-6])
    gmsh.model.geo.addPlaneSurface([cl])

    cl = gmsh.model.geo.addCurveLoop([-8,10,-3,9])
    gmsh.model.geo.addPlaneSurface([cl])

    cl = gmsh.model.geo.addCurveLoop([-16,-12,-7,11])
    gmsh.model.geo.addPlaneSurface([cl])

    cl = gmsh.model.geo.addCurveLoop([17,15,-13,14])
    gmsh.model.geo.addPlaneSurface([cl])

    cl = gmsh.model.geo.addCurveLoop([-18,19,16,-15])
    gmsh.model.geo.addPlaneSurface([cl])

    thicklines = [4,6,9,10,2,7,16,18,14,15,19,21]
    steglines = [11,12]
    flanschlines = [1,5,3,8,13,17,20,22]

    faces = [1,2,3,4,5,6,7]


    facedimtags = []
    for i in faces:
        facedimtags.append((2,i))

    outDimTags = gmsh.model.geo.extrude(facedimtags,L,0,0,[nl],recombine=True)

    print(outDimTags)
    print(outDimTags[1:len(outDimTags):6])

    frontfaceDimTags = outDimTags[0:len(outDimTags):6]
    volumeDimTags = outDimTags[1:len(outDimTags):6]
    frontfaces = []
    for i in frontfaceDimTags:
        frontfaces.append(i[1])

    volumeTags = []
    for i in volumeDimTags:
        volumeTags.append(i[1])

    gmsh.model.geo.synchronize()

    for i in thicklines:
        gmsh.model.mesh.setTransfiniteCurve(i,nthick)

    for i in steglines:
        gmsh.model.mesh.setTransfiniteCurve(i,nsteg)

    for i in flanschlines:
        gmsh.model.mesh.setTransfiniteCurve(i,nflansch)

    for i in faces:
        gmsh.model.mesh.setTransfiniteSurface(i)


    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.generate(3)

    gmsh.model.mesh.removeDuplicateNodes()
    #gmsh.model.mesh.removeDuplicateElements()

    gmsh.model.mesh.reverse(facedimtags)

    #gmsh.fltk.run()


    pathname = os.path.dirname(sys.argv[0])
    currPath = os.path.abspath(pathname)
    fesys= HierAMuS.FEMPy(pathname, "IProfCantilever")
    fesys.setStaticSolutionState()
    fesys.setSolver(4)

    mesh = fesys.getMeshCommands()
    gm = mesh.getFromGMESH()
    geo = mesh.getGeometryCommands()
    macro = fesys.getMacroCommands()

    macro.setLogLevel(fesys.FullLog(),fesys.FullLog())

    gm.addGeomFromGmsh(gmsh)
    geo.checkGeometry()

    gm.addBrickVolumeFiniteElements(gmsh,1,volumeTags,1)

    mesh.getElementFormulations().addEL300_3DSolid(1,1,disporder,1)
    mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1,100,0.3)
    mesh.addMaterial(1,1,1)


    nt,cn,t =gmsh.model.mesh.getNodes(dim=0,tag=p1)
    vnum = nt[0]
    fnums = gm.getFaceNumbers(gmsh,frontfaces,4,1)
    matlist = []
    for i in fnums:
        matlist.append(1)
    mesh.getElementCommands().add3DInterfaceElement(materialNumber=2, specialNumber=1, vertNum=vnum, faceList=fnums,faceMaterialList=matlist)
    mesh.getElementFormulations().addEL302_BeamCoupling3D(num=2, disporder=disporder, meshiddisp=1, meshidrot=2,mode=3,warptype=warptype)
    mesh.addMaterial(2, 1, 2)

    nt,cn,t =gmsh.model.mesh.getNodes(dim=0,tag=p2)
    vnum1 = nt[0]
    fnums = gm.getFaceNumbers(gmsh,faces,4,1)
    matlist = []
    for i in fnums:
        matlist.append(1)
    mesh.getElementCommands().add3DInterfaceElement(materialNumber=3, specialNumber=2, vertNum=vnum1, faceList=fnums,faceMaterialList=matlist)
    mesh.getElementFormulations().addEL302_BeamCoupling3D(num=3, disporder=disporder, meshiddisp=1, meshidrot=2,mode=3,warptype=warptype)
    mesh.addMaterial(3, 1, 3)


    mesh.setDegreesOfFreedom()

    #fnums = gm.getFaceNumbers(gmsh,faces,4,1)
    #print(faces)
    #mesh.getBoundaryConditions().singleBC(geo.faceType(),fnums,1,[1,1,1],shapeOrder=disporder)
    #fnums = gm.getFaceNumbers(gmsh,frontfaces,4,1)
    #print(frontfaces)
    #mesh.getBoundaryConditions().singleLoad(geo.faceType(),fnums,1,[0,0,1],1)

    mesh.getBoundaryConditions().singleBC(geo.vertexType(), vnum1, 1, [1,1,1], shapeOrder=disporder)
    mesh.getBoundaryConditions().singleBC(geo.vertexType(), vnum1, 2, [1,1,1], shapeOrder=disporder)
    mesh.getBoundaryConditions().singleLoad(geo.vertexType(),vnum,2,[0,1,0],1,shapeorder=disporder)

    macro.sparseSetUp()

    macro.setPropFunction(1)
    macro.setDt(1)
    macro.timeincr()

    macro.assembleSolve()
    macro.computeEigenValues(10,20,max=False)
    #macro.computeEigenValues(3,9,max=True)

    fesys.getPlotCommands().toFile()



run()