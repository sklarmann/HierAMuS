# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import gmsh
import HierAMuS
import math
import numpy

def createGMSHGeomtry():
    
    L=4
    b=1
    h=1
    ll=0.03
    
    nw = 8
    nh = 8
    nL = 10
    
    # Load vector
    loadv = numpy.array([0.0,1.,0])
    
    # Rotation axis
    ax = numpy.array([1.0,0.0,0.0])
    ax /= numpy.linalg.norm(ax)
    # Rotation angle
    theta = 1.5*math.pi/4
    theta = 0
    
    loadv = math.cos(theta)*loadv + math.sin(theta)*numpy.cross(ax,loadv) + (1-math.cos(theta))*numpy.dot(ax,loadv)*ax
    
    gmsh.initialize()
    gmsh.model.add('rectbeam')
    
    p1 = gmsh.model.occ.addPoint(-L/2,0,0)
    p2 = gmsh.model.occ.addPoint(L/2,0,0)
    
    
    
    gmsh.model.occ.addBox(-L/2+ll,-b/2,-h/2,L-2*ll,b,h)
    
    gmsh.model.occ.rotate([(0,2),(3,1)],-L/2,0,0,ax[0],ax[1],ax[2],theta)
    
    gmsh.model.occ.synchronize()
    #gmsh.fltk.run()
    
    lengthlinetags = [9,10,11,12]
    widthlinetags = [2,4,6,8]
    heightlinetags = [1,3,5,7]
    
    facetags = [1,2,3,4,5,6]
    
    
    
    gmsh.model.occ.synchronize()
    
    for i in widthlinetags:
        gmsh.model.mesh.setTransfiniteCurve(i,nw)
        
    for i in heightlinetags:
        gmsh.model.mesh.setTransfiniteCurve(i,nh)
    
    for i in lengthlinetags:
        gmsh.model.mesh.setTransfiniteCurve(i,nL)
    
    for i in facetags:
        gmsh.model.mesh.setTransfiniteSurface(i)
    
    gmsh.model.mesh.setTransfiniteVolume(1)
    
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    
    gmsh.model.mesh.generate(3)
    
    
    #gmsh.model.mesh.reverse([(2,1)])
    #gmsh.model.mesh.reverse([(2,2)])
    #gmsh.fltk.run()
    
    
    
    
    meshidDisp = 1
    meshidRot = 2
    warptype = 1
    order = 2
    interfacemode = 3
    
    
    import os, sys
    pathname = os.path.dirname(sys.argv[0])
    
    fesys = HierAMuS.FEMPy(pathname,"RectangularPlotTest")
    fesys.setStaticSolutionState()
    fesys.setSolver(4)
    
    
    macro = fesys.getMacroCommands()
    mesh = fesys.getMeshCommands()
    geo = mesh.getGeometryCommands()
    gm = mesh.getFromGMESH()
    
    macro.setLogLevel(fesys.BasicLog(),fesys.BasicLog())
    
    gm.addGeomFromGmsh(gmsh)
    geo.checkGeometry()
    
    gm.addVolumeElements(gmsh,1,1)
    
    mesh.getElementFormulations().addEL300_3DSolid(num=1,meshiddisp=meshidDisp,disporder=order,mode=1)
    mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(number=1,E=100,nu=0.3)
    mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(number=2,E=100,nu=0.0)
    mesh.addMaterial(matNum=1,matFormNum=1,elemFormNum=1)
    
    
    # left interface
    fnums = gm.getFaceNumbers(gmsh,1,4,1)
    fmat = []
    for i in fnums:
        fmat.append(1)
    mesh.getElementCommands().add3DInterfaceElement(materialNumber=2,specialNumber=1,vertNum=1,faceList=fnums,faceMaterialList=fmat)
    mesh.getElementFormulations().addEL302_BeamCoupling3D(num=2,disporder=order,meshiddisp=meshidDisp,meshidrot=meshidRot,mode=1,warptype=warptype)
    mesh.addMaterial(matNum=2,matFormNum=1,elemFormNum=2)
    
    # right interface
    fnums = gm.getFaceNumbers(gmsh,2,4,1)
    fmat = []
    for i in fnums:
        fmat.append(1)
    mesh.getElementCommands().add3DInterfaceElement(materialNumber=3,specialNumber=2,vertNum=2,faceList=fnums,faceMaterialList=fmat)
    mesh.getElementFormulations().addEL302_BeamCoupling3D(num=3,disporder=order,meshiddisp=meshidDisp,meshidrot=meshidRot,mode=4,warptype=1)#,warpBounNodes=[24,12,20])
    mesh.addMaterial(matNum=3,matFormNum=1,elemFormNum=3)
    
    mesh.setDegreesOfFreedom()
    bc = mesh.getBoundaryConditions()
    
    #fnum = gm.getFaceNumbers(gmsh,1,4,1)
    #bc.singleBC(eltype=geo.faceType(),number=fnum,meshId=meshidDisp,dofs=[1,1,1],shapeOrder=order)
    
    
    bc.BCVertex(number=2,meshId=meshidDisp,dofs=[1,1,1])
    bc.BCVertex(number=2,meshId=meshidRot,dofs=[1,1,1])
    #bc.BCVertex(number=2,meshId=meshidRot,dofs=[1,1,1])
    
    #bc.LoadVertex(number=2,meshId=meshidRot,load=[loadv[0],loadv[1],loadv[2]],propnum=1)
    bc.LoadVertex(number=1,meshId=meshidDisp,load=[0,1,0],propnum=1)
    
    
    macro.sparseSetUp()
    
    macro.setPropFunction(number=1)
    macro.setDt(1)
    macro.timeincr()
    
    macro.assembleSolve()
    #macro.assembleSolve()
    
    #macro.computeEigenValues(10,40)
    
    fesys.getPlotCommands().toFile()
    soldisp = macro.getSolution(geo.vertexType(),2,1)
    solrot = macro.getSolution(geo.vertexType(),2,2)
    print("Displacement: ",soldisp)
    print("Rotation:     ",solrot)
    
createGMSHGeomtry()
