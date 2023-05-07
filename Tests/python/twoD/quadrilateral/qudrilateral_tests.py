# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import unittest
import os, sys
import HierAMuS

def firstModel(order,meshid):
    pathname = os.path.dirname(sys.argv[0])
    currPath = os.path.abspath(pathname)


    L  = 10
    h1 = 1

    nx = 5
    ny = 2



    fesys = HierAMuS.FEMPy(currPath, 'firsttest')
    fesys.setStaticSolutionState()
    fesys.setSolver(1)
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

    fesys.getMacroCommands().printInfo()
    sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(),geomNumber=v3.geoId,meshId=meshid)
    #print(sol)
    return sol

class TestQuad(unittest.TestCase):
    def __init__(self, methodName: str = "QuadTest"):
        self.places=8
        super().__init__(methodName)

    def test_quad1(self):
        res = firstModel(1,1)
        self.assertAlmostEqual(res[0], 1.1667516758203988, places=self.places, msg="Test failed for order 1")
        self.assertAlmostEqual(res[1], 15.637013564123457, places=self.places, msg="Test failed for order 1")
    def test_quad2(self):
        res = firstModel(2,2)
        self.assertAlmostEqual(res[0], 2.9909727974564415, places=self.places, msg="Test failed for order 2")
        self.assertAlmostEqual(res[1], 39.95185373540363, places=self.places, msg="Test failed for order 2")
    #     self.assertSequenceEqual(firstModel(2,2), [2.9909727974564415, 39.95185373540363, 0.0], "Test failed for order 2")
    def test_quad3(self):
        res = firstModel(3,3)
        self.assertAlmostEqual(res[0], 2.9949503288809844, places=self.places, msg="Test failed for order 3")
        self.assertAlmostEqual(res[1], 40.14916844740804, places=self.places, msg="Test failed for order 3")
    #     self.assertSequenceEqual(firstModel(3,3), [2.9949503288809844, 40.14916844740804, 0.0], "Test failed for order 3")
    def test_quad4(self):   
        res = firstModel(4,4)
        self.assertAlmostEqual(res[0], 2.9977137053220146, places=self.places, msg="Test failed for order 4")
        self.assertAlmostEqual(res[1], 40.196681228814285, places=self.places, msg="Test failed for order 4")
    #     self.assertSequenceEqual(firstModel(4,4), [2.9977137053220146, 40.196681228814285, 0.0], "Test failed for order 4")
    def test_quad5(self):
        res = firstModel(5,5)
        self.assertAlmostEqual(res[0], 2.99922379128834, places=self.places, msg="Test failed for order 5")
        self.assertAlmostEqual(res[1], 40.21849279832855, places=self.places, msg="Test failed for order 5")
    #     self.assertSequenceEqual(firstModel(5,5), [2.99922379128834, 40.21849279832855, 0.0], "Test failed for order 5")
    def test_quad6(self):
        res = firstModel(6,6)
        self.assertAlmostEqual(res[0], 3.000077517578813, places=self.places, msg="Test failed for order 6")
        self.assertAlmostEqual(res[1], 40.22914352116219, places=self.places, msg="Test failed for order 6")
    def test_quad7(self):
        res = firstModel(7,7)
        self.assertAlmostEqual(res[0], 3.0005532862061686, places=self.places, msg="Test failed for order 7")
        self.assertAlmostEqual(res[1], 40.2347592648595, places=self.places, msg="Test failed for order 7")
    def test_quad8(self):
        res = firstModel(8,8)
        self.assertAlmostEqual(res[0], 3.0008311466987783, places=self.places, msg="Test failed for order 8")
        self.assertAlmostEqual(res[1], 40.23804062904741, places=self.places, msg="Test failed for order 8")
    def test_quad9(self):
        res = firstModel(9,9)
        self.assertAlmostEqual(res[0], 3.00100810206853, places=self.places, msg="Test failed for order 9")
        self.assertAlmostEqual(res[1], 40.240127227211566, places=self.places, msg="Test failed for order 9")
    def test_quad10(self):
        res = firstModel(10,10)
        self.assertAlmostEqual(res[0], 3.001128688515539, places=self.places, msg="Test failed for order 10")
        self.assertAlmostEqual(res[1], 40.24152427921814, places=self.places, msg="Test failed for order 10")
    #     self.assertAlmostEqual(firstModel(6,6), [3.000077517578813, 40.22914352116219, 0.0], "Test failed for order 6")


if __name__ == '__main__':
    unittest.main()