# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import os, sys
import HierAMuS
import gmsh
import unittest

def firstModel(order,meshid,solver):
    pathname = os.path.dirname(sys.argv[0])
    currPath = os.path.abspath(pathname)
    
    L = 10
    h = 1
    b = 1
    h1 = 1
    nx = 2
    ny = 3
    nz = 3

    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)

    gmsh.model.add('CantileverBeam')
    
    p1=gmsh.model.occ.addPoint(0,0,0)
    p2=gmsh.model.occ.addPoint(L,0,0)
    p3=gmsh.model.occ.addPoint(L,h,0)
    p4=gmsh.model.occ.addPoint(0,h,0)
    
    p5=gmsh.model.occ.addPoint(0,0,b)
    p6=gmsh.model.occ.addPoint(L,0,b)
    p7=gmsh.model.occ.addPoint(L,h,b)
    p8=gmsh.model.occ.addPoint(0,h,b)
    
    l1 = gmsh.model.occ.addLine(p1,p2)
    l2 = gmsh.model.occ.addLine(p2,p3)
    l3 = gmsh.model.occ.addLine(p4,p3)
    l4 = gmsh.model.occ.addLine(p1,p4)
    
    lengthLines = [l1,l3]
    widthLines = [l2,l4]
    
    l5 = gmsh.model.occ.addLine(p1,p5)
    l6 = gmsh.model.occ.addLine(p2,p6)
    l7 = gmsh.model.occ.addLine(p3,p7)
    l8 = gmsh.model.occ.addLine(p4,p8)
    
    heightLines=[l5,l6,l7,l8]
    
    l9 = gmsh.model.occ.addLine(p5,p6)
    l10 = gmsh.model.occ.addLine(p6,p7)
    l11 = gmsh.model.occ.addLine(p8,p7)
    l12 = gmsh.model.occ.addLine(p5,p8)
    
    lengthLines += [l9,l11]
    widthLines += [l10,l12]
    
    
    cl1 = gmsh.model.occ.addCurveLoop([l1,l2,l3,l4])
    cl2 = gmsh.model.occ.addCurveLoop([l9,l10,l11,l12])
    
    cl3 = gmsh.model.occ.addCurveLoop([l2,l7,l10,l6])
    cl4 = gmsh.model.occ.addCurveLoop([l4,l8,l12,l5])
    
    
    cl5 = gmsh.model.occ.addCurveLoop([l1,l5,l9,l6])
    cl6 = gmsh.model.occ.addCurveLoop([l3,l8,l11,l7])
    
    f1 = gmsh.model.occ.addPlaneSurface([cl1])
    f2 = gmsh.model.occ.addPlaneSurface([cl2])
    
    f3 = gmsh.model.occ.addPlaneSurface([cl3])
    f4 = gmsh.model.occ.addPlaneSurface([cl4])
    
    f5 = gmsh.model.occ.addPlaneSurface([cl5])
    f6 = gmsh.model.occ.addPlaneSurface([cl6])
    
    faces = [f1,f2,f3,f4,f5,f6]
    
    fl1 = gmsh.model.occ.addSurfaceLoop([f1,f2,f3,f4,f5,f6])
    v1 = gmsh.model.occ.addVolume([fl1])
    gmsh.model.occ.synchronize()
    
    for i in lengthLines:
        gmsh.model.mesh.setTransfiniteCurve(i,nx+1)
    
    for i in widthLines:
        gmsh.model.mesh.setTransfiniteCurve(i,ny+1)
    
    for i in heightLines:
        gmsh.model.mesh.setTransfiniteCurve(i,nz+1)
    
    for i in faces:
        gmsh.model.mesh.setTransfiniteSurface(i)
        
    gmsh.model.mesh.setTransfiniteVolume(v1)
    
    
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.model.mesh.generate(3)
    
    gmsh.fltk.run()


    fesys = HierAMuS.FEMPy(currPath, 'firsttest')
    fesys.getMacroCommands().setLogLevel(fesys.NoLog(),fesys.NoLog())
    fesys.setStaticSolutionState()
    fesys.setSolver(solver)
    
    geo = fesys.getMeshCommands().getGeometryCommands()
    
    fromGM = fesys.getMeshCommands().getFromGMESH()
    fromGM.addGeomFromGmsh(gmsh)
    geo.checkGeometry()
    
    # Adding Element
    fromGM.addVolumeElements(gmsh,v1,1)
    
    # Element definition
    fesys.getMeshCommands().getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1,100,0.3)
    fesys.getMeshCommands().getElementFormulations().addEL300_3DSolid(1,1,order,1)
    fesys.getMeshCommands().addMaterial(1,matFormNum=1,elemFormNum=1)

    fesys.getMeshCommands().setDegreesOfFreedom()
    
    bfaces = fromGM.getFaceNumbers(gmsh,f3,4,1)
    lfaces = fromGM.getFaceNumbers(gmsh,f4,4,1)
    fesys.getMeshCommands().getBoundaryConditions().BC(geo.faceType(),bfaces,1,[1,1,1],order)
    
    
    fesys.getMeshCommands().getBoundaryConditions().load(geo.faceType(),lfaces,1,[0,1,0],0)
    
    print(lfaces)
    print(bfaces)
    gmsh.fltk.run()
    
    fesys.getMacroCommands().sparseSetUp()

    fesys.getMacroCommands().setPropFunction(number=0)
    fesys.getMacroCommands().setDt(1)
    fesys.getMacroCommands().timeincr()

    fesys.getMacroCommands().newton(refResidual=1e-11)
    
    sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(),geomNumber=p6,meshId=meshid)

    fesys.getPlotCommands().toFile()
    
    return sol

firstModel(1,1,1)
# class TestConstraint2D(unittest.TestCase):
#     def __init__(self, methodName: str = "Constraint2D"):
#         self.places=8
#         super().__init__(methodName)
        
#     def test_const1(self):
#         sol=firstModel(1,1,1)
#         self.assertAlmostEqual(sol[0],1.8756191968421505,places=self.places, msg="Test failed for solver 1")
#         self.assertAlmostEqual(sol[1],12.733402020529594,places=self.places, msg="Test failed for solver 1")
#     def test_const2(self):
#         sol=firstModel(1,1,2)
#         self.assertAlmostEqual(sol[0],1.8756191968421505,places=self.places, msg="Test failed for solver 2")
#         self.assertAlmostEqual(sol[1],12.733402020529594,places=self.places, msg="Test failed for solver 2")
#     def test_const3(self):
#         sol=firstModel(1,1,3)
#         self.assertAlmostEqual(sol[0],1.8756191968421505,places=self.places, msg="Test failed for solver 3")
#         self.assertAlmostEqual(sol[1],12.733402020529594,places=self.places, msg="Test failed for solver 3")
#     def test_const4(self):
#         sol=firstModel(1,1,4)
#         self.assertAlmostEqual(sol[0],1.8756191968421505,places=self.places, msg="Test failed for solver 4")
#         self.assertAlmostEqual(sol[1],12.733402020529594,places=self.places, msg="Test failed for solver 4")
#     def test_const5(self):
#         sol=firstModel(1,1,5)
#         self.assertAlmostEqual(sol[0],1.8756191968421505,places=self.places, msg="Test failed for solver 5")
#         self.assertAlmostEqual(sol[1],12.733402020529594,places=self.places, msg="Test failed for solver 5")
#     def test_const6(self):
#         sol=firstModel(1,1,6)
#         self.assertAlmostEqual(sol[0],1.8756191968421505,places=self.places, msg="Test failed for solver 6")
#         self.assertAlmostEqual(sol[1],12.733402020529594,places=self.places, msg="Test failed for solver 6")

