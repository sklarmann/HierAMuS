# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
import gmsh
import HierAMuS
import os,sys
import math


def run(disporder,geoorder,nr,nphi):
    
    
    
    rl = 2
    b=4
            
    gmsh.initialize()
    gmsh.model.add('lochscheibe')
    gmsh.option.setNumber('General.Terminal', 0)
    
    
    cc=gmsh.model.occ.addPoint(0,0,0)
    
    p1=gmsh.model.occ.addPoint(rl,0,0)
    p2=gmsh.model.occ.addPoint(b,0,0)
    p3=gmsh.model.occ.addPoint(b,b,0)
    p4=gmsh.model.occ.addPoint(rl*math.cos(45/180*math.pi),rl*math.sin(45/180*math.pi),0)
    p5=gmsh.model.occ.addPoint(0,rl,0)
    p6=gmsh.model.occ.addPoint(0,b,0)
    
    circleLine1=gmsh.model.occ.addCircleArc(p4,cc,p1)
    circleLine2=gmsh.model.occ.addCircleArc(p4,cc,p5)
    bottomLine = gmsh.model.occ.addLine(p1,p2)
    l4 = gmsh.model.occ.addLine(p2,p3)
    centerLine = gmsh.model.occ.addLine(p3,p4)
    topLine = gmsh.model.occ.addLine(p3,p6)
    leftLine = gmsh.model.occ.addLine(p6,p5)
    
    cl1 = gmsh.model.occ.addCurveLoop([circleLine1,bottomLine,l4,centerLine])
    f1 = gmsh.model.occ.addPlaneSurface([cl1])
    cl2 = gmsh.model.occ.addCurveLoop([circleLine2,centerLine,topLine,leftLine])
    f2 = gmsh.model.occ.addPlaneSurface([cl2])
    
    
    gmsh.model.occ.synchronize()
    
    radialLines = [bottomLine,centerLine,leftLine]
    otherLines = [circleLine1,l4,topLine,circleLine2]
    
    for i in radialLines:
        gmsh.model.mesh.setTransfiniteCurve(i,nphi)
        
    for i in otherLines:
        gmsh.model.mesh.setTransfiniteCurve(i,nr)
    
    gmsh.model.mesh.setTransfiniteSurface(f1)
    gmsh.model.mesh.setTransfiniteSurface(f2)
    
    gmsh.option.setNumber("Mesh.Algorithm", 4)
    gmsh.option.setNumber("Mesh.RecombineAll", 2)
    
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(geoorder)
    gmsh.model.mesh.reverse([(2,f2)])
    
    #gmsh.fltk.run()
    
    
    path = os.path.dirname(sys.argv[0])
    print(__file__)
    fesys = HierAMuS.FEMPy(path,"scheibeloch")
    fesys.setStaticSolutionState()
    fesys.setSolver(3)
    fesys.getMacroCommands().setLogLevel(fesys.NoLog(),fesys.NoLog())
    
    gm = fesys.getMeshCommands().getFromGMESH()
    geo = fesys.getMeshCommands().getGeometryCommands()
    gm.addGeomFromGmsh(gmsh)
    geo.checkGeometry()
    
    gm.addFaceElements(gmsh,f1,1)
    gm.addFaceElements(gmsh,f2,1)
    fesys.getMeshCommands().getElementFormulations().addEL201_2DShell(num=1,meshiddisp=1,disporder=disporder,mode=1)
    fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(number=1,E=2100,nu=0.3,thickness=1,plainstrain=0)
    fesys.getMeshCommands().addMaterial(1,1,1)
    
    fesys.getMeshCommands().setDegreesOfFreedom()
    lnums = gm.getEdgeNumbers(gmsh,bottomLine,geoorder)
    fesys.getMeshCommands().getBoundaryConditions().singleBC(eltype=geo.edgeType(),number=lnums,meshId=1,dofs=[0,1,1],shapeOrder=disporder)
    lnums = gm.getEdgeNumbers(gmsh,leftLine,geoorder)
    fesys.getMeshCommands().getBoundaryConditions().singleBC(eltype=geo.edgeType(),number=lnums,meshId=1,dofs=[1,0,1],shapeOrder=disporder)
    
    lnums = gm.getEdgeNumbers(gmsh,topLine,geoorder)
    #print(lnums)
    fesys.getMeshCommands().getBoundaryConditions().singleLoad(eltype=geo.edgeType(),number=lnums,meshId=1,load=[0,1,0],propnum=1,shapeorder=disporder)
    
    fesys.getMacroCommands().sparseSetUp()
    
    fesys.getMacroCommands().setPropFunction(1)
    fesys.getMacroCommands().setDt(1)
    fesys.getMacroCommands().timeincr()
    
    fesys.getMacroCommands().assembleSolve()
    
    [coor,parmcoor,dim,nnode] = gmsh.model.mesh.getNode(p3)
    
    sol = fesys.getMacroCommands().getSolution(geomType=geo.vertexType(),geomNumber=nnode,meshId=1)
    #fesys.getPlotCommands().toFile()
    gmsh.finalize()
    return sol
    
    
nr = 2
nphi = 2
nn=5

x=[]
y1=[]
y2=[]

for i in range(nn):
    
    sol = run(3,1,nr,nphi)
    y1.append(sol[1])
    sol = run(3,2,nr,nphi)
    y2.append(sol[1])
    x.append(i)
    print(sol)
    nr*=2
    nphi*=2
    
print(x,y1)
    
fig, ax = plt.subplots()

ax.plot(x, y1, linewidth=3.0)
ax.plot(x, y2, linewidth=1.0)
#ax.plot(numref, quadsol, linewidth=2.0)

plt.show()