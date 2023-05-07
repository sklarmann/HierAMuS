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
    h1 = 1
    nx = 2
    ny = 3

    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)

    gmsh.model.add('CantileverBeam')
    
    p1=gmsh.model.occ.addPoint(0,0,0)
    p2=gmsh.model.occ.addPoint(L/2,0,0)
    p3=gmsh.model.occ.addPoint(L/2,h,0)
    p4=gmsh.model.occ.addPoint(0,h,0)
    
    
    p5=gmsh.model.occ.addPoint(L/2,0,0)
    p6=gmsh.model.occ.addPoint(L,0,0)
    p7=gmsh.model.occ.addPoint(L,h,0)
    p8=gmsh.model.occ.addPoint(L/2,h,0)
    
    l1=gmsh.model.occ.addLine(p1,p2)
    l2=gmsh.model.occ.addLine(p2,p3)
    l3=gmsh.model.occ.addLine(p3,p4)
    l4=gmsh.model.occ.addLine(p4,p1)
    
    l5=gmsh.model.occ.addLine(p5,p6)
    l6=gmsh.model.occ.addLine(p6,p7)
    l7=gmsh.model.occ.addLine(p7,p8)
    l8=gmsh.model.occ.addLine(p8,p5)
    
    cl1=gmsh.model.occ.addCurveLoop([l1,l2,l3,l4])
    f1=gmsh.model.occ.addPlaneSurface([cl1])
    cl2=gmsh.model.occ.addCurveLoop([l5,l6,l7,l8])
    f2=gmsh.model.occ.addPlaneSurface([cl2])

    gmsh.model.occ.synchronize()
    
    gmsh.model.mesh.setTransfiniteCurve(l1,nx+1)
    gmsh.model.mesh.setTransfiniteCurve(l3,nx+1)
    gmsh.model.mesh.setTransfiniteCurve(l2,ny+1)
    gmsh.model.mesh.setTransfiniteCurve(l4,ny+1)
    gmsh.model.mesh.setTransfiniteSurface(f1)
    
    gmsh.model.mesh.setTransfiniteCurve(l5,nx+1)
    gmsh.model.mesh.setTransfiniteCurve(l7,nx+1)
    gmsh.model.mesh.setTransfiniteCurve(l6,ny+1)
    gmsh.model.mesh.setTransfiniteCurve(l8,ny+1)
    gmsh.model.mesh.setTransfiniteSurface(f2)
    
    
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.model.mesh.generate(4)
    


    fesys = HierAMuS.FEMPy(currPath, 'firsttest')
    fesys.getMacroCommands().setLogLevel(fesys.NoLog(),fesys.NoLog())
    fesys.setStaticSolutionState()
    fesys.setSolver(solver)
    
    geo = fesys.getMeshCommands().getGeometryCommands()
    
    fromGM = fesys.getMeshCommands().getFromGMESH()
    fromGM.addGeomFromGmsh(gmsh)
    geo.checkGeometry()
    
    # Adding Element
    fromGM.addFaceElements(gmsh,[f1,f2],1)
    
    # Element definition
    fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(1,E=100,nu=0.3,thickness=1,plainstrain=0)
    fesys.getMeshCommands().getElementFormulations().addEL201_2DShell(num=1,meshiddisp=meshid,disporder=order,mode=1)
    fesys.getMeshCommands().addMaterial(1,matFormNum=1,elemFormNum=1)

    fesys.getMeshCommands().setDegreesOfFreedom()
    
    
    
    # Boundary Conditions
    edgeList = fromGM.getEdgeNumbers(gmsh,l4,1)
    fesys.getMeshCommands().getBoundaryConditions().BC(eltype=geo.edgeType(),number=edgeList,meshId=1,dofs=[1,1,1],shapeOrder=order)

    masterEdges = fromGM.getEdgeNumbers(gmsh,l2,1)
    slaveEdges = fromGM.getEdgeNumbers(gmsh,l8,1)
    a=[1,2,3]
    a.reverse()
    slaveEdges.reverse()
    

    fesys.getMeshCommands().getConstraintCommands().generalLink(geo.edgeType(),masterEdges,slaveEdges,1,order,0,0,1,1)
    fesys.getMeshCommands().getConstraintCommands().generalLink(geo.edgeType(),masterEdges,slaveEdges,1,order,1,1,1.0,1.0)

    loadEdges = fromGM.getEdgeNumbers(gmsh,l6,1)
    fesys.getMeshCommands().getBoundaryConditions().load(geo.edgeType(),loadEdges,1,[0,1,0],0)

    fesys.getMacroCommands().sparseSetUp()

    fesys.getMacroCommands().setPropFunction(number=0)
    fesys.getMacroCommands().setDt(1)
    fesys.getMacroCommands().timeincr()

    fesys.getMacroCommands().newton(refResidual=1e-11)
    
    sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(),geomNumber=p6,meshId=meshid)

    
    return sol


class TestConstraint2D(unittest.TestCase):
    def __init__(self, methodName: str = "Constraint2D"):
        self.places=8
        super().__init__(methodName)
        
    def test_const1(self):
        sol=firstModel(1,1,1)
        self.assertAlmostEqual(sol[0],1.8756191968421505,places=self.places, msg="Test failed for solver 1")
        self.assertAlmostEqual(sol[1],12.733402020529594,places=self.places, msg="Test failed for solver 1")
    def test_const2(self):
        sol=firstModel(1,1,2)
        self.assertAlmostEqual(sol[0],1.8756191968421505,places=self.places, msg="Test failed for solver 2")
        self.assertAlmostEqual(sol[1],12.733402020529594,places=self.places, msg="Test failed for solver 2")
    def test_const3(self):
        sol=firstModel(1,1,3)
        self.assertAlmostEqual(sol[0],1.8756191968421505,places=self.places, msg="Test failed for solver 3")
        self.assertAlmostEqual(sol[1],12.733402020529594,places=self.places, msg="Test failed for solver 3")
    def test_const4(self):
        sol=firstModel(1,1,4)
        self.assertAlmostEqual(sol[0],1.8756191968421505,places=self.places, msg="Test failed for solver 4")
        self.assertAlmostEqual(sol[1],12.733402020529594,places=self.places, msg="Test failed for solver 4")
    def test_const5(self):
        sol=firstModel(1,1,5)
        self.assertAlmostEqual(sol[0],1.8756191968421505,places=self.places, msg="Test failed for solver 5")
        self.assertAlmostEqual(sol[1],12.733402020529594,places=self.places, msg="Test failed for solver 5")
    def test_const6(self):
        sol=firstModel(1,1,6)
        self.assertAlmostEqual(sol[0],1.8756191968421505,places=self.places, msg="Test failed for solver 6")
        self.assertAlmostEqual(sol[1],12.733402020529594,places=self.places, msg="Test failed for solver 6")

