# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import os, sys
import HierAMuS

def firstModel(order,meshid):
    pathname = os.path.dirname(sys.argv[0])
    currPath = os.path.abspath(pathname)


    L  = 4
    h1 = 4

    nx = 1
    ny = 1



    fesys = HierAMuS.FEMPy(currPath, 'firsttest')
    fesys.setStaticHomogenizationSolutionState()
    fesys.setSolver(2)
    fesys.getMacroCommands().setLogLevel(fesys.BasicLog(),fesys.BasicLog())
    geoB = fesys.getGeomBuilder()

    # Lower Block
    # Vertices definieren
    v1 = geoB.getVertex()
    v1.setCoordinates(-L/2,-h1/2,0)
    v2 = geoB.getVertex()
    v2.setCoordinates(L/2,-h1/2,0)
    v3 = geoB.getVertex()
    v3.setCoordinates(L/2,h1/2,0)
    v4 = geoB.getVertex()
    v4.setCoordinates(-L/2,h1/2,0)
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
    
    fesys.getMeshCommands().getGeometryCommands().checkGeometry()

    fesys.getMeshCommands().getElementCommands().addFace(1,q1.getQuadList())

    fesys.getMeshCommands().getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1,E=100,nu=0.0)
    fesys.getMeshCommands().getElementFormulations().addEL300_3DSolid(num=1,meshiddisp=meshid,disporder=order,mode=1)
    #fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(1,E=100,nu=0.0,thickness=1,plainstrain=0)
    #fesys.getMeshCommands().getElementFormulations().addEL201_2DShell(num=1,meshiddisp=meshid,disporder=order,mode=1)
    
    fesys.getMeshCommands().addMaterial(1,matFormNum=1,elemFormNum=1)

    fesys.getMeshCommands().setDegreesOfFreedom()
    
    fesys.getMacroCommands().getHomogenizationCommands().setHomogenizationSolid3D(meshIdDisp=1,dispOrder=order,bctype=0)
    
    fesys.getMacroCommands().sparseSetUp()
    
    fesys.getMacroCommands().getHomogenizationCommands().computeAMatrix()

    #fesys.getMacroCommands().setPropFunction(number=0)
    #fesys.getMacroCommands().setDt(1)
    #fesys.getMacroCommands().timeincr()
    
    fesys.getMacroCommands().getHomogenizationCommands().setStrains([0.0,0.02,0.0])

    fesys.getMacroCommands().newton()
    #fesys.getMacroCommands().assembleSolve()
    #fesys.getMacroCommands().assembleSolve()
    fesys.getMacroCommands().getHomogenizationCommands().homogenize()

    sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(),geomNumber=v3.geoId,meshId=meshid)
    print(sol)
    fesys.getPlotCommands().toFile()
    #fesys.getMacroCommands().setLogLevel(fesys.FullLog(),fesys.FullLog())
    #fesys.getMacroCommands().printInfo()
    return sol

firstModel(1,1)