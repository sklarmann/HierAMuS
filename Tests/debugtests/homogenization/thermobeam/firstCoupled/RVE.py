# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import os,sys
import HierAMuS
import gmsh


def rectAngularRVE(L,b,h,nx,ny,nz,order):
    # current path of file
    currPath = os.path.dirname(os.path.realpath(__file__))


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
    RVE.setSolver(4)
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
    
    mform = RVE.getMeshCommands().getMaterialFormulations()
    #mform.addMA1_3D_LinearElastic_Isotrop(1,100,0.3)
    mform.addMA3_SmallStrainPlasticity(number=1,E=100,nu=0.3,y0=10,yinf=0,xh=10,xd=0,eta=0)
    
    eform = RVE.getMeshCommands().getElementFormulations()
    eform.addEL300_3DSolid(num=1,meshiddisp=1,disporder=order,mode=1)
    eform.addEL307_VolumeConstraint(num=2,meshiddisp=1,meshIdLam=3,shapeorder=order,center=1,stiffness=1)
    eform.addEL207_FaceConstraint(num=3,meshIdDisp=1,meshIdRot=5,meshIdLam=6,meshIdMu=7,dispOrder=order,k=1,mode=1)
    
    mesh = RVE.getMeshCommands()
    mesh.addMaterial(1,1,1)
    mesh.addMaterial(2,1,2)
    mesh.addMaterial(3,1,3)
    
    mesh.setDegreesOfFreedom()
    mesh.getBoundaryConditions().BCVertex(cVert,1,[1,1,1])
    mesh.getBoundaryConditions().BCVertex(cVert,5,[1,1,1])
    
    RVE.getMacroCommands().getHomogenizationCommands().setHomogenizationBeam(1,order,1)
    
    RVE.getMacroCommands().sparseSetUp()
    
    
    RVE.getMacroCommands().getHomogenizationCommands().computeAMatrix()
    
    gmsh.finalize()
    
    return RVE


def run():
    #pass
    return 0
    RVE = rectAngularRVE(2,1,1,4,4,4,order=2)
    RVE.getMacroCommands().getHomogenizationCommands().setStrains([0.1,0.0,0.0,0.0,0.0,0.0])
    RVE.getMacroCommands().assembleSolve()
    RVE.getMacroCommands().assembleSolve()
    RVE.getMacroCommands().getHomogenizationCommands().homogenize()
    
    
run()