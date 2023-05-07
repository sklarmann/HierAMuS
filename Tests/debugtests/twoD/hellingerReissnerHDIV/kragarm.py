# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
import gmsh
import HierAMuS
import os, sys


def run(refinements,geometryOrder,disporder,quad=False,basic=False):
    gmsh.initialize()
    gmsh.model.add('test')
    gmsh.option.setNumber('General.Terminal', 0)
    
    L=10
    h=1
    p1=gmsh.model.occ.addPoint(0,0,0)
    p2=gmsh.model.occ.addPoint(L,0,0)
    p3=gmsh.model.occ.addPoint(L,h,0)
    p4=gmsh.model.occ.addPoint(0,h,0)

    l1=gmsh.model.occ.addLine(p1,p2)
    l2=gmsh.model.occ.addLine(p2,p3)
    l3=gmsh.model.occ.addLine(p3,p4)
    l4=gmsh.model.occ.addLine(p4,p1)
    
    cl = gmsh.model.occ.addCurveLoop([l1,l2,l3,l4])
    f1 = gmsh.model.occ.addPlaneSurface([cl])
    
    gmsh.model.occ.synchronize()
    
    gmsh.option.setNumber("Mesh.Algorithm", 5)
    if quad is True:
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
    else:
        gmsh.option.setNumber("Mesh.RecombineAll", 0)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 0)
    gmsh.option.setNumber("Mesh.ElementOrder", geometryOrder)
    
    gmsh.model.mesh.generate(2)
    
    for i in range(refinements):
        gmsh.model.mesh.refine()
    
    #gmsh.fltk.run()
    
    path = os.path.dirname(sys.argv[0])
    
    fesys = HierAMuS.FEMPy(path,"cantilever")
    fesys.setStaticSolutionState()
    fesys.setSolver(4)
    fesys.getMacroCommands().setLogLevel(fesys.NoLog(),fesys.NoLog())
    
    gm = fesys.getMeshCommands().getFromGMESH()
    geo = fesys.getMeshCommands().getGeometryCommands()
    gm.addGeomFromGmsh(gmsh)
    geo.checkGeometry()
    
    gm.addFaceElements(gmsh,f1,1)
    if basic is True:
        fesys.getMeshCommands().getElementFormulations().addEL201_2DShell(num=1,meshiddisp=1,disporder=disporder,mode=1)
    else:
        fesys.getMeshCommands().getElementFormulations().addEL205_HDivTest(num=1,plainstrain=0,disporder=disporder,stressorder=disporder-1,mode=2,meshiddisp=1,meshidstress=2,E=100,nu=0.3)
    fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(number=1,E=100,nu=0.3,thickness=1,plainstrain=0)
    fesys.getMeshCommands().addMaterial(1,1,1)
    
    fesys.getMeshCommands().setDegreesOfFreedom()
    lnums = gm.getEdgeNumbers(gmsh,l4,geometryOrder)
    
    fesys.getMeshCommands().getBoundaryConditions().singleBC(eltype=geo.edgeType(),number=lnums,meshId=1,dofs=[1,1,1],shapeOrder=disporder)
    
    lnums = gm.getEdgeNumbers(gmsh,l2,geometryOrder)

    fesys.getMeshCommands().getBoundaryConditions().singleLoad(eltype=geo.edgeType(),number=lnums,meshId=1,load=[0,1,0],propnum=1,shapeorder=disporder)
    
    fesys.getMacroCommands().sparseSetUp()
    
    fesys.getMacroCommands().setPropFunction(1)
    fesys.getMacroCommands().setDt(1)
    fesys.getMacroCommands().timeincr()
    
    fesys.getMacroCommands().assembleSolve()
    
    [coor,parmcoor,dim,nnode] = gmsh.model.mesh.getNode(p3)
    
    #fesys.getPlotCommands().toFile()
    return fesys.getMacroCommands().getSolution(geomType=geo.vertexType(),geomNumber=nnode,meshId=1)
    


x=[]

basicQuadOrder1 = []
basicQuadOrder2 = []
basicTriOrder1 = []
basicTriOrder2 = []
HRQuadOrder1 = []
HRQuadOrder2 = []

nn=7
for i in range(nn):
    x.append(i)
    sol = run(refinements=i+1,geometryOrder=1,disporder=1,quad=False,basic=True)
    basicTriOrder1.append(sol[1])
    print(sol)
    
    sol = run(refinements=i,geometryOrder=1,disporder=2,quad=False,basic=True)
    basicTriOrder2.append(sol[1])
    print(sol)
    
    sol = run(refinements=i+1,geometryOrder=1,disporder=1,quad=True,basic=True)
    basicQuadOrder1.append(sol[1])
    print(sol)
    
    sol = run(refinements=i,geometryOrder=1,disporder=2,quad=True,basic=True)
    basicQuadOrder2.append(sol[1])
    print(sol)
    
    sol = run(refinements=i+1,geometryOrder=1,disporder=1,quad=True,basic=False)
    HRQuadOrder1.append(sol[1])
    print(sol)
    
    sol = run(refinements=i,geometryOrder=1,disporder=2,quad=True,basic=False)
    HRQuadOrder2.append(sol[1])
    print(sol)
    
    
fig, ax = plt.subplots()
#fig, ax2 = plt.subplots()


ll, =ax.plot(x, basicTriOrder1, linewidth=1.0)
ll.set_label('TriBasic 1')

ll, =ax.plot(x, basicTriOrder2, linewidth=1.0)
ll.set_label('TriBasic 2')


ll, =ax.plot(x, basicQuadOrder1, linewidth=1.0)
ll.set_label('QuadBasic 1')

ll, =ax.plot(x, basicQuadOrder2, linewidth=1.0)
ll.set_label('QuadBasic 2')

ll, =ax.plot(x, HRQuadOrder1, linewidth=1.0)
ll.set_label('HRQuadOrder 1')

ll, =ax.plot(x, HRQuadOrder2, linewidth=1.0)
ll.set_label('HRQuadOrder 2')
#ax.plot(numref, quadsol, linewidth=2.0)

ax.legend()
ax.set_ylim(39.8,40.3)
#ax2.legend()
plt.show()