# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

# from xml.sax.handler import feature_string_interning
import os
import sys

import HierAMuS
import gmsh

def run(quadMesh,order,numrefs):

    gmsh.initialize()
    gmsh.model.add('test')


    L=10
    h=1

    E=100
    nu=0.0
    G=E/2/(1+nu)


    disporder=order

    p1=gmsh.model.occ.addPoint(0,0,0,2)
    p2=gmsh.model.occ.addPoint(L,0,0,2)
    p3=gmsh.model.occ.addPoint(L,h,0,2)
    p4=gmsh.model.occ.addPoint(0,h,0,2)

    l1=gmsh.model.occ.addLine(p1,p2)
    l2=gmsh.model.occ.addLine(p2,p3)
    l3=gmsh.model.occ.addLine(p3,p4)
    l4=gmsh.model.occ.addLine(p4,p1)

    cl=gmsh.model.occ.addCurveLoop([l1,l2,l3,l4])
    face=gmsh.model.occ.addPlaneSurface([cl])

    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.Algorithm",4)

    if quadMesh:
        gmsh.option().setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)

    gmsh.model.mesh.generate(2)
    for i in range(numrefs):
        gmsh.model.mesh.refine()

    gmsh.fltk.run()

    pathname = os.path.dirname(sys.argv[0])
    currPath = os.path.abspath(pathname)
    fesys= HierAMuS.FEMPy(pathname, "firstvolume")
    fesys.setStaticSolutionState()
    fesys.setSolver(3)
    fesys.getMacroCommands().setLogLevel(fesys.NoLog(), fesys.NoLog())

    gm = fesys.getMeshCommands().getFromGMESH()
    gm.addGeomFromGmsh(gmsh)

    fesys.getMeshCommands().getGeometryCommands().checkGeometry()

    gm.addFaceElements(gmsh,face,1)

    # if quadMesh:
    #     gm.addQuadrilateralFiniteElements(gmsh,1,face,1)
    # else:
    #     gm.addTriangleFiniteElements(gmsh,1,face,1)


    fesys.getMeshCommands().getElementFormulations().addEL201_2DShell(1,meshiddisp=1,disporder=disporder,mode=1)
    #fesys.getMeshCommands().getElementFormulations().addEL205_HDivTest(1,100,disporder=disporder,stressorder=disporder-1,mode=2,meshiddisp=1,meshidstress=2,E=E,nu=nu)
    fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(1,E=E,nu=nu,thickness=1,plainstrain=0)
    fesys.getMeshCommands().addMaterial(1,1,1)

    fesys.getMeshCommands().setDegreesOfFreedom()

    edges = gm.getEdgeNumbers(gmsh,l4,1)
    fesys.getMeshCommands().getBoundaryConditions().BC(fesys.getMeshCommands().getGeometryCommands().edgeType(),edges,meshId=1,dofs=[1,1,1],shapeOrder=disporder)
    #fesys.getMeshCommands().getBoundaryConditions().BC(fesys.getMeshCommands().getGeometryCommands().vertexType(),p1,meshId=1,dofs=[1,1,1],shapeOrder=disporder)
    #fesys.getMeshCommands().getBoundaryConditions().BC(fesys.getMeshCommands().getGeometryCommands().vertexType(),p4,meshId=1,dofs=[1,0,1],shapeOrder=disporder)
    
    edges = gm.getEdgeNumbers(gmsh,l2,1)
    fesys.getMeshCommands().getBoundaryConditions().load(fesys.getMeshCommands().getGeometryCommands().edgeType(),edges,meshId=1,load=[0,1,0],propnum=0,shapeorder=disporder)

    fesys.getMacroCommands().sparseSetUp()

    fesys.getMacroCommands().setPropFunction(0)
    fesys.getMacroCommands().setDt(1)
    fesys.getMacroCommands().timeincr()

    fesys.getMacroCommands().assembleSolve()

    fesys.getPlotCommands().toFile()


    # [nt,ncoor,npcoor]=gmsh.model.mesh.getNodes(1,l2,includeBoundary=True)

    # t=0
    # for i in nt:
    #     sol=fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(),i,meshId=1)
    #     print(sol)
    #     t+=sol[1]

    # print(t/len(nt))

    # print("Result should be: ",4*L*L*L/E/h/h/h+L/G*6/5)
    
    sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(),p2,1)
    fesys.getMacroCommands().setLogLevel(fesys.FullLog(), fesys.BasicLog())
    fesys.getMacroCommands().printInfo()
    fesys.getMacroCommands().computeEigenValues(10,30,max=False)
    return sol


print(run(False,2,1))