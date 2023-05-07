# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import os,sys
import HierAMuS
import gmsh
import matplotlib.pyplot as plt
import csv

def run(L,nx,ny,nz,order,bctype):
    gmsh.initialize()
    gmsh.model.add("test")

    currPath = os.path.dirname(os.path.realpath(__file__))
    geoFile = os.path.join(currPath,"RVE.geo")

    # FE1
    E=100
    nu=0.3
    alpha = 1
    #alpha = 22.8*10**-6
    #alpha=0
    #c = 452
    c = 0
    rho0 = 0
    T0 = 0
    kappa = 100

    lam = E*nu/((1+nu)*(1-2*nu))
    G=E/(2*(1+nu))


    b=1
    h=1



    #gmsh.open(geoFile)




    gmsh.model.occ.addPoint(-L/2, -b/2, -h/2, 1.0)

    gmsh.model.occ.addPoint(-L/2, b/2, -h/2, 1.0)
    gmsh.model.occ.addPoint(-L/2, b/2, h/2, 1.0)
    gmsh.model.occ.addPoint(-L/2, -b/2, h/2, 1.0)
    gmsh.model.occ.addPoint(0, -b/2, -h/2, 1.0)
    gmsh.model.occ.addPoint(0, b/2, -h/2, 1.0)
    gmsh.model.occ.addPoint(0, b/2, h/2, 1.0)
    gmsh.model.occ.addPoint(0, -b/2, h/2, 1.0)
    gmsh.model.occ.addPoint(L/2, -b/2, -h/2, 1.0)
    gmsh.model.occ.addPoint(L/2, b/2, -h/2, 1.0)
    gmsh.model.occ.addPoint(L/2, b/2, h/2, 1.0)
    gmsh.model.occ.addPoint(L/2, -b/2, h/2, 1.0)
    gmsh.model.occ.addPoint(0, 0, 0, 1.0)
    gmsh.model.occ.addPoint(0, 0, 0, 1.0)

    gmsh.model.occ.addLine(4, 8)
    gmsh.model.occ.addLine(8, 5)
    gmsh.model.occ.addLine(5, 1)
    gmsh.model.occ.addLine(1, 4)
    gmsh.model.occ.addLine(4, 3)
    gmsh.model.occ.addLine(3, 2)
    gmsh.model.occ.addLine(2, 1)
    gmsh.model.occ.addLine(2, 6)
    gmsh.model.occ.addLine(6, 10)
    gmsh.model.occ.addLine(10, 11)
    gmsh.model.occ.addLine(11, 7)
    gmsh.model.occ.addLine(7, 6)
    gmsh.model.occ.addLine(7, 3)
    gmsh.model.occ.addLine(7, 8)
    gmsh.model.occ.addLine(5, 6)
    gmsh.model.occ.addLine(10, 9)
    gmsh.model.occ.addLine(9, 5)
    gmsh.model.occ.addLine(9, 12)
    gmsh.model.occ.addLine(12, 11)
    gmsh.model.occ.addLine(12, 8)

    gmsh.model.occ.addCurveLoop([5, 6, 7, 4])
    gmsh.model.occ.addPlaneSurface([1])
    gmsh.model.occ.addCurveLoop([8, -15, 3, -7])
    gmsh.model.occ.addPlaneSurface([2])
    gmsh.model.occ.addCurveLoop([12, -8, -6, -13])
    gmsh.model.occ.addPlaneSurface([3])
    gmsh.model.occ.addCurveLoop([3, 4, 1, 2])
    gmsh.model.occ.addPlaneSurface([4])
    gmsh.model.occ.addCurveLoop([1, -14, 13, -5])
    gmsh.model.occ.addPlaneSurface([5])
    gmsh.model.occ.addCurveLoop([14, 2, 15, -12])
    gmsh.model.occ.addPlaneSurface([6])
    gmsh.model.occ.addCurveLoop([20, 2, -17, 18])
    gmsh.model.occ.addPlaneSurface([7])
    gmsh.model.occ.addCurveLoop([17, 15, 9, 16])
    gmsh.model.occ.addPlaneSurface([8])
    gmsh.model.occ.addCurveLoop([12, 9, 10, 11])
    gmsh.model.occ.addPlaneSurface([9])
    gmsh.model.occ.addCurveLoop([11, 14, -20, 19])
    gmsh.model.occ.addPlaneSurface([10])
    gmsh.model.occ.addCurveLoop([18, 19, -10, 16])
    gmsh.model.occ.addPlaneSurface([11])

    gmsh.model.occ.synchronize()

    for i in [3, 1, 13, 8, 11, 20, 9, 17]:
        gmsh.model.mesh.setTransfiniteCurve(i,nx)
    for i in [4, 6, 2, 12, 18, 10]:
        gmsh.model.mesh.setTransfiniteCurve(i,nz)
    for i in [5, 7, 14, 15, 19, 16]:
        gmsh.model.mesh.setTransfiniteCurve(i,ny)



    gmsh.model.mesh.setTransfiniteSurface(1)

    gmsh.model.mesh.setTransfiniteSurface(4)

    gmsh.model.mesh.setTransfiniteSurface(5)

    gmsh.model.mesh.setTransfiniteSurface(3)

    gmsh.model.mesh.setTransfiniteSurface(2)

    gmsh.model.mesh.setTransfiniteSurface(6)

    gmsh.model.mesh.setTransfiniteSurface(10)

    gmsh.model.mesh.setTransfiniteSurface(9)

    gmsh.model.mesh.setTransfiniteSurface(8)

    gmsh.model.mesh.setTransfiniteSurface(7)

    gmsh.model.mesh.setTransfiniteSurface(11)



    gmsh.model.occ.addSurfaceLoop([5, 4, 2, 3, 1, 6])

    gmsh.model.occ.addVolume([1])

    gmsh.model.occ.addSurfaceLoop([10, 9, 8, 7, 11, 6])

    gmsh.model.occ.addVolume([2])

    gmsh.model.occ.synchronize()

    gmsh.model.mesh.setTransfiniteVolume(1)
    gmsh.model.mesh.setTransfiniteVolume(2)

    gmsh.model.addPhysicalGroup(3,[1,2],tag=21,name="volelements")
    gmsh.model.addPhysicalGroup(2,[6],22,name="interface")

    gmsh.model.addPhysicalGroup(0,[13],23,"interfacePoint")
    gmsh.model.addPhysicalGroup(0,[14],24,"interfacePoint")


    gmsh.option.setNumber("Mesh.RecombineAll", 1) # This option sets gmsh to recombine tetra to bricks

    gmsh.model.mesh.generate(3)

    #gmsh.fltk.run()

    Rve = HierAMuS.FEMPy(currPath,"ThermoRVEPeriodic")
    Rve.setStaticHomogenizationSolutionState()
    Rve.setSolver(6)
    Rve.getMacroCommands().setLogLevel(Rve.NoLog(),Rve.NoLog())

    mesh = Rve.getMeshCommands()
    macro = Rve.getMacroCommands()
    gm = mesh.getFromGMESH()
    geo = mesh.getGeometryCommands()


    gm.addGeomFromGmsh(gmsh)
    geo.checkGeometry()



    v1 = gmsh.model.getEntitiesForPhysicalGroup(3,21)
    #gmsh.fltk.run()



    gm.addVolumeElements(gmsh,v1.tolist(),1)
    mesh.getElementFormulations().addEL303_ThermoMechanikSolid3D(num=1,meshiddisp=1,meshidtemperatur=2,shapeorder=order,mu=G,lamb=lam,alpha=alpha,c=c,rho0=rho0,T0=T0,kappa=kappa,mode=1)
    #mesh.getElementFormulations().addEL300_3DSolid(num=1,meshiddisp=1,disporder=order,mode=1)
    mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(number=1,E=E,nu=nu)
    mesh.addMaterial(matNum=1,matFormNum=1,elemFormNum=1)


    cent = 0
    ft = gmsh.model.getEntitiesForPhysicalGroup(2,22).tolist()
    nt = gmsh.model.getEntitiesForPhysicalGroup(0,23).tolist()[0]
    if bctype==1:
        gm.addFaceConstraint(gmsh,faceTags=ft,vertexTag=nt,material=2)
        mesh.getElementFormulations().addEL207_FaceConstraint(2,1,1,3,4,5,dispOrder=order,k=1)
        mesh.addMaterial(2,1,2)
        cent=1

    
    nt2 = gmsh.model.getEntitiesForPhysicalGroup(0,24).tolist()[0]
    gm.addVolumeConstraint(gmsh,v1.tolist(),nt,3)
    mesh.getElementFormulations().addEL307_VolumeConstraint(num=3,meshiddisp=1,meshIdLam=10,shapeorder=order,stiffness=1,center=cent)
    mesh.addMaterial(3,1,3)


    #nt2 = gmsh.model.getEntitiesForPhysicalGroup(0,24).tolist()[0]
    #gm.addVolumeConstraint(gmsh,v1.tolist(),nt,4)
    #mesh.getElementFormulations().addEL307_VolumeConstraint(num=4,meshiddisp=2,meshIdLam=11,shapeorder=order,stiffness=1,center=0,mode=4)
    #mesh.addMaterial(4,1,4)


    mesh.setDegreesOfFreedom()
    if bctype == 1:
        mesh.getBoundaryConditions().BC(geo.vertexType(),nt,3,[1,1,1],1)
        mesh.getBoundaryConditions().BC(geo.vertexType(),nt,1,[1,1,1],1)

    macro.getHomogenizationCommands().Homogenization3DThermoMechBeam(meshIdDisp=1,dispOrder=order,meshIdTemp=2,tempOrder=order,bctype=bctype)

    macro.sparseSetUp()

    macro.getHomogenizationCommands().computeAMatrix()


    eps = 0
    deps = 1
    macro.setDt(1)

    eps=0
    for i in range(1):
        macro.timeincr()
        eps+=deps
        macro.getHomogenizationCommands().setStrains([0,0,0,0,0,0,0,0,0,1,0,0])


        macro.newton(refResidual=1e-9)
        #Rve.getPlotCommands().toFile()

    macro.getHomogenizationCommands().homogenize()

    Rve.getPlotCommands().toFile()
    return macro.getHomogenizationCommands().getCMatrix()



vals = {"EA" : [], "GA" : [], "GI" : [], "EI" : [], "EAalph" : [], "EIalph" : [], "kappaA" : [], "kappaEI" : [], "x":[]}


steps=1
Lmod=1
x=[]
for i in range(steps):
    Lmod*=2
    x.append(Lmod)
    C = run(L=Lmod,nx=40,ny=7,nz=7,order=2,bctype=1)
    vals["x"].append(Lmod)
    vals["EA"].append(C[0][0])
    vals["GA"].append(C[1][1])
    vals["GI"].append(C[3][3])
    vals["EI"].append(C[4][4])
    vals["EAalph"].append(C[0][6])
    vals["EIalph"].append(C[4][7])
    vals["kappaA"].append(C[7][7])
    vals["kappaEI"].append(C[10][10])
    

#Lmod = 0.01
#vals2 = {"EA" : [], "GA" : [], "GI" : [], "EI" : [], "EAalph" : [], "EIalph" : [], "kappaA" : [], "kappaEI" : [], "x": []}
#for i in range(steps):
#    Lmod*=2
#    x.append(Lmod)
#    C = run(L=Lmod,nx=4,ny=6,nz=6,order=2,bctype=0)
#    vals2["x"].append(Lmod)
#    vals2["EA"].append(C[0][0])
#    vals2["GA"].append(C[1][1])
#    vals2["GI"].append(C[3][3])
#    vals2["EI"].append(C[4][4])
#    vals2["EAalph"].append(C[0][6])
#    vals2["EIalph"].append(C[4][7])
#    vals2["kappaA"].append(C[7][7])
#    vals2["kappaEI"].append(C[10][10])

#print(vals)


fig, ax = plt.subplots(1,1)
#ax.plot(vals["x"],vals["EI"])
#print(vals["x"])
#print(vals["EA"])




currPath = os.path.dirname(os.path.realpath(__file__))
csvFile = os.path.join(currPath,"data.csv")

with open(csvFile,'w') as csvfile:
    for nn, var in vals.items():
        csvfile.write(nn)
        for j in var:
            csvfile.write(",")
            csvfile.write(str(j))
        csvfile.write("\n")
        
        

#plt.show()