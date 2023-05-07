# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import os, sys
import HierAMuS

def firstModel(order,meshid):
    pathname = os.path.dirname(sys.argv[0])
    currPath = os.path.abspath(pathname)


    L  = 10
    h1 = 1

    nx = 2
    ny = 2



    fesys = HierAMuS.FEMPy(currPath, 'firsttest')
    fesys.setStaticHomogenizationSolutionState()
    fesys.setSolver(2)
    geoB = fesys.getGeomBuilder()

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

    fesys.getMeshCommands().getBoundaryConditions().singleBC(fesys.getMeshCommands().getGeometryCommands().edgeType(),number=e4.EdgeList,meshId=meshid,dofs=[1,1,1],shapeOrder=order)
    fesys.getMeshCommands().getBoundaryConditions().singleLoad(fesys.getMeshCommands().getGeometryCommands().edgeType(),e2.EdgeList,meshId=meshid,load=[0,1,0],propnum=0,shapeorder=order)


    fesys.getMacroCommands().sparseSetUp()

    fesys.getMacroCommands().setPropFunction(number=0)
    fesys.getMacroCommands().setDt(1)
    fesys.getMacroCommands().timeincr()

    fesys.getMacroCommands().assembleSolve()

    sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(),geomNumber=v3.geoId,meshId=meshid)
    print(sol)
    fesys.getMacroCommands().setLogLevel(fesys.FullLog(),fesys.FullLog())
    
    fesys.getMacroCommands().solutionStateToFile("solutiontest.bin")
    fesys.getMacroCommands().solutionStateFromFile("solutiontest.bin")
    fesys.getMacroCommands().timeincr()
    fesys.getMacroCommands().assembleSolve()

    fesys.getPlotCommands().toFile()
    #fesys.getMacroCommands().printInfo()
    return sol

firstModel(1,1)