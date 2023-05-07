# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import os, sys
import HierAMuS
import gmsh



def run():

    E=210000*10**6 # N/m^2
    v=0.3

    #E=1
    #v = 0

    lamb = E*v/(1+v)/(1-2*v)
    mu = E/2/(1+v)

    alpha = 22.8*10**-6 # 1/K
    #alpha = 0
    kappa = 81 # W/m/K
    kappa = 0
    cv = 452 # J/kg/K
    #cv = 0
    rho0 = 7860 # kg/m^3

    #print(lamb, mu, alpha, kappa)
    #exit()
    gmsh.initialize()

    nlT = 3
    ncT = 3
    nb  = 3
    nl  = 6

    ll = 0.0
    L = 10
    L-=2*ll

    disporder = 2

    gmsh.model.add("RectangularBeam")

    gmsh.model.occ.addPoint(ll, -0.5, -0.5)
    gmsh.model.occ.addPoint(ll,  0.5, -0.5)
    gmsh.model.occ.addPoint(ll,  0.5, 0.5)
    gmsh.model.occ.addPoint(ll, -0.5, 0.5)



    gmsh.model.occ.addLine(1, 2)
    gmsh.model.occ.addLine(2, 3)
    gmsh.model.occ.addLine(3, 4)
    gmsh.model.occ.addLine(4, 1)



    lines = [1,2,3,4]
    gmsh.model.occ.addCurveLoop(lines)

    gmsh.model.occ.addPlaneSurface([1])

    gmsh.model.occ.synchronize()
    gmsh.fltk.run()


    backfaces = [1]


    gmsh.model.occ.extrude([(2,1)], L, 0, 0,[nl],recombine=True)
    gmsh.model.occ.synchronize()
    gmsh.fltk.run()


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


    gmsh.model.occ.synchronize()

    gmsh.fltk().run()


    pathname = os.path.dirname(sys.argv[0])
    currPath = os.path.abspath(pathname)
    fesys= HierAMuS.FEMPy(pathname, "RectangularCantilever")
    fesys.setStaticSolutionState()
    fesys.setSolver(3)
    

    mesh = fesys.getMeshCommands()
    gm = mesh.getFromGMESH()
    geo = mesh.getGeometryCommands()
    macro = fesys.getMacroCommands()

    macro.setLogLevel(fesys.FullLog(),fesys.FullLog())


    gm.addGeomFromGmsh(gmsh)
    geo.checkGeometry()

    gm.addBrickVolumeFiniteElements(gmsh, 1, [1], 1)

    mesh.getElementFormulations().addEL303_ThermoMechanikSolid3D(1,1,2,disporder,mu=mu,lamb=lamb,alpha=alpha,c=cv,rho0=rho0,T0=300,kappa=kappa,mode=1)
    mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1,100,0.0)
    mesh.addMaterial(1,1,1)

    mesh.setDegreesOfFreedom()

    fnums = gm.getFaceNumbers(gmsh,1,4,1)
    mesh.getBoundaryConditions().BC(geo.faceType(), fnums, 1, [1,1,1], shapeOrder=disporder)
    #mesh.getBoundaryConditions().singleBC(geo.faceType(), fnums, 2, [1,1,1], shapeOrder=disporder)

    #fnums = gm.getFaceNumbers(gmsh,5,4,1)
    #[nt, temp1, temp2] = gmsh.model.mesh.getNodes(dim=2,tag = 5,includeBoundary=True)
    #mesh.getBoundaryConditions().singleBC(geo.faceType(), fnums, 2, [1,1,1], shapeOrder=disporder)
    #mesh.getBoundaryConditions().singleLoad(geo.vertexType(), nt.tolist(), 2, [290,0,0], 1)
#
    #fnums = gm.getFaceNumbers(gmsh,3,4,1)
    #[nt, temp1, temp2] = gmsh.model.mesh.getNodes(dim=2,tag = 3,includeBoundary=True)
    #mesh.getBoundaryConditions().singleBC(geo.faceType(), fnums, 2, [1,1,1], shapeOrder=disporder)
    #mesh.getBoundaryConditions().singleLoad(geo.vertexType(), nt.tolist(), 2, [310,0,0], 1)


    fnums = gm.getFaceNumbers(gmsh,6,4,1)
    mesh.getBoundaryConditions().load(geo.faceType(), fnums, 1, [0,100,0], 1)

    #fnums = gm.getFaceNumbers(gmsh,backfaces,4,1)
    #mesh.getBoundaryConditions().singleLoad(geo.faceType(), fnums, 1, [0,1,0], 1)

    macro.sparseSetUp()

    macro.setPropFunction(1)
    macro.setDt(1)
    macro.timeincr()

    macro.newton()
    #macro.computeEigenValues(10,40)
    #macro.computeEigenValues(10,40,max=True)

    fesys.getPlotCommands().toFile()




run()