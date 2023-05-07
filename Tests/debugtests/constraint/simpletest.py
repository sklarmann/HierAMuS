# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import os, sys
import HierAMuS


def solvePrint(fesys,v2,v3,meshid):
    fesys.getMacroCommands().assembleSolve()
    print("Printing solution")
    sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(),geomNumber=v3.num,meshId=meshid)
    print(sol)
    sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(),geomNumber=v2.num,meshId=meshid)
    print(sol)

def firstModel(order,meshid):
    pathname = os.path.dirname(sys.argv[0])
    currPath = os.path.abspath(pathname)


    L  = 10
    h1 = 1

    nx = 4
    ny = 4



    fesys = HierAMuS.FEMPy(currPath, 'firsttest')
    fesys.setStaticSolutionState()
    fesys.setSolver(2)
    geoB = fesys.getGeomBuilder()
    geo = fesys.getMeshCommands().getGeometryCommands()
    fesys.getMacroCommands().setLogLevel(fesys.BasicLog(),fesys.BasicLog())

    # Lower Block
    # Vertices definieren
    v1 = geoB.getVertex()
    v1.setCoordinates(0,0,0)
    v2 = geoB.getVertex()
    v2.setCoordinates(L,0,0)
    v3 = geoB.getVertex()
    v3.setCoordinates(L,h1,0)
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

    fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(1,E=100,nu=0.3,thickness=1,plainstrain=0)
    fesys.getMeshCommands().getElementFormulations().addEL201_2DShell(num=1,meshiddisp=meshid,disporder=order,mode=1)
    fesys.getMeshCommands().addMaterial(1,matFormNum=1,elemFormNum=1)

    fesys.getMeshCommands().setDegreesOfFreedom()
    fesys.getMeshCommands().getBoundaryConditions().BC(eltype=geo.edgeType(),number=e4.EdgeList,meshId=1,dofs=[1,1,1],shapeOrder=order)

    fesys.getMeshCommands().getConstraintCommands().generalLink(geoType=geo.vertexType(),masterNumber=v3.num,slaveNumber=v2.num,meshId=1,masterDof=1,slaveDof=1,factor=1,difference=0.4)
    fesys.getMeshCommands().getConstraintCommands().generalLink(geoType=geo.vertexType(),masterNumber=v3.num,slaveNumber=v2.num,meshId=1,masterDof=0,slaveDof=0,factor=1,difference=0.2)
    
    #fesys.getMeshCommands().getBoundaryConditions().load(eltype=geo.vertexType(),number=v3.num,meshId=1,load=[0,1,0],propnum=0)

    fesys.getMacroCommands().sparseSetUp()

    fesys.getMacroCommands().setPropFunction(number=0)
    fesys.getMacroCommands().setDt(1)
    fesys.getMacroCommands().timeincr()

    fesys.getMacroCommands().newton(refResidual=1e-11)
    #fesys.getMacroCommands().assembleSolve()
    #fesys.getMacroCommands().assembleSolve()
    print(" ")
    solvePrint(fesys,v2,v3,meshid)
    solvePrint(fesys,v2,v3,meshid)
    print(" ")

    #sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(),geomNumber=v3.num,meshId=meshid)
    #print(sol)
    sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(),geomNumber=v2.num,meshId=meshid)
    #print(sol)
    
    fesys.getMacroCommands().solutionStateToFile("solutiontest.bin")
    fesys.getMacroCommands().solutionStateFromFile("solutiontest.bin")
    #fesys.getMacroCommands().timeincr()
    #fesys.getMacroCommands().newton(refResidual=1e-11)
    #fesys.getMacroCommands().assembleSolve()
    #fesys.getMacroCommands().assembleSolve()

    fesys.getPlotCommands().toFile()
    #fesys.getMacroCommands().printInfo()
    print("calculation done")
    return sol

firstModel(1,1)