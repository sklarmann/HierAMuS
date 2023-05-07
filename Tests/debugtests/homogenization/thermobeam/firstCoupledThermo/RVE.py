# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import os,sys
import HierAMuS
import gmsh
import matplotlib.pyplot as plt


def rectAngularRVE(L,b,h,nx,ny,nz,order):
    # current path of file
    currPath = os.path.dirname(os.path.realpath(__file__))

    E=100
    nu=0.3
    
    
    lam = E*nu/((1+nu)*(1-2*nu))
    mu = E/(2*(1+nu))

    gmsh.initialize()
    gmsh.model.add("RVE")


    cVert = gmsh.model.occ.addPoint(0,0,0)
    gmsh.model.occ.addBox(-L/2, -b/2, -h/2, L/2, b, h, 1)
    gmsh.model.occ.addBox(0, -b/2, -h/2, L/2,b, h, 2)

    
    gmsh.model.occ.synchronize()
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    lengthlines = [9,10,11,12,17,18,19,20]
    heightlines = [1,3,5,7,13,15]
    widhtlines = [2,4,6,8,14,16]

    for i in lengthlines:
        gmsh.model.mesh.setTransfiniteCurve(i,nx+1)
    for j in heightlines:
        gmsh.model.mesh.setTransfiniteCurve(j,nz+1)
    for k in widhtlines:
        gmsh.model.mesh.setTransfiniteCurve(k,ny+1)

    faces = [1,2,3,4,5,6,7,8,9,10,11]
    for i in faces:
        gmsh.model.mesh.setTransfiniteSurface(i)    

    cFace = 2

    gmsh.model.mesh.setTransfiniteVolume(1)
    gmsh.model.mesh.setTransfiniteVolume(2)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.model.mesh.generate(3)
    

    RVE = HierAMuS.FEMPy(currPath,"RVE")
    RVE.setStaticHomogenizationSolutionState()
    RVE.setSolver(6)
    RVE.getMacroCommands().setLogLevel(RVE.BasicLog(),RVE.BasicLog())
    
    gm = RVE.getMeshCommands().getFromGMESH()
    gm.addGeomFromGmsh(gmsh)
    #gmsh.fltk.run()
    
    RVE.getMeshCommands().getGeometryCommands().checkGeometry()
    
    gm.addVolumeElements(gmsh,1,1)
    gm.addVolumeElements(gmsh,2,1)
    gm.addVolumeConstraint(gmsh,1,cVert,2)
    gm.addVolumeConstraint(gmsh,2,cVert,2)
    gm.addFaceConstraint(gmsh,cFace,cVert,3)
    #gm.addVolumeConstraint(gmsh,1,cVert,4)
    #gm.addVolumeConstraint(gmsh,2,cVert,4)
    
    mform = RVE.getMeshCommands().getMaterialFormulations()
    mform.addMA1_3D_LinearElastic_Isotrop(1,100,0.3)
    #mform.addMA3_SmallStrainPlasticity(number=1,E=100,nu=0.3,y0=10,yinf=0,xh=10,xd=0,eta=0)
    
    eform = RVE.getMeshCommands().getElementFormulations()
    #eform.addEL300_3DSolid(num=1,meshiddisp=1,disporder=order,mode=1)
    eform.addEL303_ThermoMechanikSolid3D(num=1,meshiddisp=1,meshidtemperatur=2,shapeorder=order,mu=mu,lamb=lam,alpha=0.1,c=0,rho0=0,T0=0,kappa=1,mode=1)
    eform.addEL307_VolumeConstraint(num=2,meshiddisp=1,meshIdLam=3,shapeorder=order,center=1,stiffness=1)
    eform.addEL207_FaceConstraint(num=3,meshIdDisp=1,meshIdRot=4,meshIdLam=5,meshIdMu=6,dispOrder=order,k=1,mode=1)
    eform.addEL307_VolumeConstraint(num=4,meshiddisp=2,meshIdLam=7,shapeorder=order,center=0,stiffness=1,mode=4)
    
    
    mesh = RVE.getMeshCommands()
    mesh.addMaterial(1,1,1)
    mesh.addMaterial(2,1,2)
    mesh.addMaterial(3,1,3)
    mesh.addMaterial(4,1,4)
    
    mesh.setDegreesOfFreedom()
    mesh.getBoundaryConditions().BCVertex(cVert,1,[1,1,1])
    mesh.getBoundaryConditions().BCVertex(cVert,4,[1,1,1])
    
    RVE.getMacroCommands().getHomogenizationCommands().Homogenization3DThermoMechBeam(meshIdDisp=1,dispOrder=order,meshIdTemp=2,tempOrder=order,bctype=1)
    
    RVE.getMacroCommands().sparseSetUp()
    
    
    RVE.getMacroCommands().getHomogenizationCommands().computeAMatrix()
    #gmsh.fltk.run()
    gmsh.finalize()
    
    return RVE


def run():
    #pass
    #return 0
    L=1.0
    RVE = rectAngularRVE(1,1,1,8,8,8,order=1)
    RVE.getMacroCommands().setLogLevel(RVE.FullLog(),RVE.FullLog())
    #RVE.getMacroCommands().printInfo()
    RVE.getMacroCommands().getHomogenizationCommands().setStrains([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.1])
    #RVE.getMacroCommands().newton(refResidual=1e-7)
    RVE.getMacroCommands().assembleSolve()
    RVE.getMacroCommands().assembleSolve()
    RVE.getMacroCommands().getHomogenizationCommands().homogenize()
    print("test")
    C=RVE.getMacroCommands().getHomogenizationCommands().getCMatrix()
    print(C[7,7])
    print(C[11,11])
    RVE.getPlotCommands().toFile()
    RVE.getMacroCommands().computeEigenValues(10,30)
    
    
    
def createDiag():
    L=0.01
    x = []
    y = []
    for i in range(16):
        RVE = rectAngularRVE(L,1,1,4,4,4,order=2)
        RVE.getMacroCommands().setLogLevel(RVE.NoLog(),RVE.NoLog())
        #RVE.getMacroCommands().printInfo()
        RVE.getMacroCommands().getHomogenizationCommands().setStrains([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.1])
        #RVE.getMacroCommands().newton(refResidual=1e-7)
        RVE.getMacroCommands().assembleSolve()
        RVE.getMacroCommands().assembleSolve()
        RVE.getMacroCommands().getHomogenizationCommands().homogenize()
        C=RVE.getMacroCommands().getHomogenizationCommands().getCMatrix()
        x.append(L)
        y.append(abs(C[11,11]))
        #RVE.getPlotCommands().toFile()
        #RVE.getMacroCommands().computeEigenValues(10,30)
        L*=2
        
    fig, ax = plt.subplots()
    ax.plot(x,y)
    plt.show()
    
    
run()