# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import os
import sys

import HierAMuS
import gmsh

quadMesh = True
HDiv = True

refinements = 5
order = 2



if HDiv == True:
    ElFormulation = "HDiv"
else:
    ElFormulation = "H1"
if quadMesh == True:
    ElType = "Q"
else:
    ElType = "T"

numEl = 4
L=48
h1=44
h2=60
meshid=1
load = 10
Emodul = 200
nu=0.499

for i in range(refinements):
    pathname = os.path.dirname(sys.argv[0])
    currPath = os.path.abspath(pathname)
    fesys = HierAMuS.FEMPy(currPath + '\log_files', f"Cook_{ElFormulation}_{ElType}_{numEl}_{order-1}")
    esize=1/numEl*L
    gmsh.initialize()
    gmsh.model.add("t11")
    p1 = gmsh.model.geo.addPoint(0,0,0,esize)
    p2 = gmsh.model.geo.addPoint(L,h1,0,esize)
    p3 = gmsh.model.geo.addPoint(L,h2,0,esize)
    p4 = gmsh.model.geo.addPoint(0,h1,0,esize)
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    face = gmsh.model.geo.addPlaneSurface([cl])

    gmsh.model.geo.synchronize()

    
    gmsh.model.mesh.setTransfiniteCurve(l1,numEl+1)
    gmsh.model.mesh.setTransfiniteCurve(l2,numEl+1)
    gmsh.model.mesh.setTransfiniteCurve(l3,numEl+1)
    gmsh.model.mesh.setTransfiniteCurve(l4,numEl+1)

    gmsh.model.mesh.setTransfiniteSurface(face, "Right")
    if quadMesh:
        gmsh.option.setNumber("Mesh.Algorithm", 8)
        gmsh.option().setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)

    gmsh.model.mesh.generate(2)

    # for j in range(i):
    #     gmsh.model.mesh.refine()

    gmsh.fltk.run()

    fesys.getMacroCommands().setLogLevel(fesys.BasicLog(), fesys.BasicLog())
    fesys.setStaticSolutionState()

    gm = fesys.getMeshCommands().getFromGMESH()
    gm.addGeomFromGmsh(gmsh)

    fesys.getMeshCommands().getGeometryCommands().checkGeometry()

    if quadMesh:
        gm.addQuadrilateralFiniteElements(gmsh,1,face,1)
    else:
        gm.addTriangleFiniteElements(gmsh,1,face,1)


    fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(1, E=Emodul, nu=nu, thickness=1, plainstrain=1)
    if HDiv:
        fesys.getMeshCommands().getElementFormulations().addEL205_HDivTest(1,1,disporder=order,stressorder=order-1,mode=2,meshiddisp=meshid,meshidstress=meshid+1,E=Emodul,nu=nu)
    else:
        fesys.getMeshCommands().getElementFormulations().addEL201_2DShell(1, meshiddisp=1, disporder=order, mode=1)


    fesys.getMeshCommands().addMaterial(1,1,1)

    fesys.getMeshCommands().setDegreesOfFreedom()


    bedges = gm.getEdgeNumbers(gmsh, l4)
    geo = fesys.getMeshCommands().getGeometryCommands()
    fesys.getMeshCommands().getBoundaryConditions().singleBC(geo.edgeType(), bedges, 1, [1,1,1], order)
    bedges = gm.getEdgeNumbers(gmsh, l2)
    fesys.getMeshCommands().getBoundaryConditions().singleLoad(geo.edgeType(),bedges,1,[0,load,0],0,shapeorder=order)


    fesys.getMacroCommands().sparseSetUp()
    fesys.getMacroCommands().setPropFunction(0)
    fesys.getMacroCommands().setDt(1)
    fesys.getMacroCommands().timeincr()
    fesys.getMacroCommands().assembleSolve()
    # fesys.getPlotCommands().toFile()

    # L2StressNorm = fesys.getMacroCommands().computeNorm(1)
    # print('Schubdiff = ' + str(L2StressNorm))

    numEls = fesys.getNumberOfElements()
    sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(),3,meshId=meshid)
    completeName = os.path.join(currPath, "Ergebnisse_Cooktest.txt")     
    f = open(completeName, "a")
    f.write(f"{f'{ElType}_numEl_{numEl}' : <12}{f'order_{order}_{ElFormulation}' : ^15}{f'{sol[0]},' : <25}{f'{sol[1]}' : <30}{f'{numEls}' : <10}\n")
    f.close()
    # gmsh.fltk.run()
    gmsh.finalize()
    #return sol
    numEl *= 2

f = open(completeName, "a")
f.write(f"\n")
f.close()