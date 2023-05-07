# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import HierAMuS
import gmsh

def run(order):
    gmsh.initialize()
    gmsh.model.add("test")
    gmsh.option.setNumber('General.Terminal', 0)

    L=10
    h=1
    b=1

    gmsh.model.occ.addPoint(0,0,0)
    gmsh.model.occ.addPoint(L,0,0)
    printNode=gmsh.model.occ.addPoint(L,b,0)
    gmsh.model.occ.addPoint(0,b,0)
    gmsh.model.occ.addPoint(0,0,h)
    gmsh.model.occ.addPoint(L,0,h)
    gmsh.model.occ.addPoint(L,b,h)
    gmsh.model.occ.addPoint(0,b,h)

    # Lines of lower face 1
    gmsh.model.occ.addLine(1,2)
    gmsh.model.occ.addLine(2,3)
    gmsh.model.occ.addLine(3,4)
    gmsh.model.occ.addLine(4,1)

    lowerFaceLines = [1,4,3,2]

    # left face liness 2
    gmsh.model.occ.addLine(1,5)
    gmsh.model.occ.addLine(5,8)
    gmsh.model.occ.addLine(8,4)

    leftFaceLines = [5,6,7,4]

    # right face lines 3
    gmsh.model.occ.addLine(2,6)
    gmsh.model.occ.addLine(6,7)
    gmsh.model.occ.addLine(7,3)

    rightFaceLines = [2,10,9,8]

    # top face lines 4
    gmsh.model.occ.addLine(5,6)
    gmsh.model.occ.addLine(7,8)

    topFaceLines = [12,6,11,9]

    # front face lines 5
    frontFaceLines = [1,8,11,5]

    # back face lines 6
    backFaceLines = [3,10,12,7]

    xlines = [1,3,11,12]
    ylines = [2,9,4,6]
    zlines = [5,7,8,10]


    # lower face 1
    cl=gmsh.model.occ.addCurveLoop(lowerFaceLines)
    f1=gmsh.model.occ.addPlaneSurface([cl])

    # left face 2
    cl=gmsh.model.occ.addCurveLoop(leftFaceLines)
    bounFace=gmsh.model.occ.addPlaneSurface([cl])


    # right face 3
    cl=gmsh.model.occ.addCurveLoop(rightFaceLines)
    loadFace=gmsh.model.occ.addPlaneSurface([cl])
    # top face 4
    cl=gmsh.model.occ.addCurveLoop(topFaceLines)
    gmsh.model.occ.addPlaneSurface([cl])
    # front face 5
    cl=gmsh.model.occ.addCurveLoop(frontFaceLines)
    gmsh.model.occ.addPlaneSurface([cl])
    # top face 6
    cl=gmsh.model.occ.addCurveLoop(backFaceLines)
    gmsh.model.occ.addPlaneSurface([cl])

    gmsh.model.occ.addSurfaceLoop([1,2,3,4,5,6])
    gmsh.model.occ.addVolume([1])


    gmsh.model.occ.synchronize()

    if order ==1:
        nx=51
        ny=5
        nz=5
    else:
        nx=26
        ny=3
        nz=3

    # settings the lines to be transfinite
    for i in xlines:
        gmsh.model.mesh.setTransfiniteCurve(i,nx)
    for i in ylines:
        gmsh.model.mesh.setTransfiniteCurve(i,ny)
    for i in zlines:
        gmsh.model.mesh.setTransfiniteCurve(i,nz)

    # setting the faces to be transfinite
    for i in [1,2,3,4,5,6]:
        gmsh.model.mesh.setTransfiniteSurface(i)

    # setting the volume to be transfinite
    gmsh.model.mesh.setTransfiniteVolume(1)


    gmsh.option.setNumber("Mesh.RecombineAll", 1) # This option sets gmsh to recombine tetra to bricks

    #gmsh.fltk.run()
    gmsh.model.mesh.generate(3)



    fesys = HierAMuS.FEMPy("./","PlasticCantilever")
    fesys.setStaticSolutionState()
    fesys.setSolver(2)

    mesh = fesys.getMeshCommands()
    macro = fesys.getMacroCommands()
    gm = mesh.getFromGMESH()
    geo = mesh.getGeometryCommands()

    macro.setLogLevel(fesys.NoLog(),fesys.NoLog())

    gm.addGeomFromGmsh(gmsh)
    geo.checkGeometry()

    gm.addVolumeElements(gmsh,1,1)

    mesh.getElementFormulations().addEL300_3DSolid(num=1,meshiddisp=1,disporder=order,mode=1)
    mesh.getMaterialFormulations().addMA3_SmallStrainPlasticity(1,E=100,nu=0.3,y0=10,yinf=0,xh=40,xd=0,eta=0)
    mesh.addMaterial(matNum=1,matFormNum=1,elemFormNum=1)

    mesh.setDegreesOfFreedom()

    fnums = gm.getFaceNumbers(gmsh=gmsh,tagin=bounFace,ftype=4,order=1)
    mesh.getBoundaryConditions().singleBC(eltype=geo.faceType(),number=fnums,meshId=1,dofs=[1,1,1],shapeOrder=order,set=True)

    fnums = gm.getFaceNumbers(gmsh=gmsh,tagin=loadFace,ftype=4,order=1)
    mesh.getBoundaryConditions().singleLoad(eltype=geo.faceType(),number=fnums,meshId=1,load=[0,1,0],propnum=1,shapeorder=order)

    [coor, pcoor, dim, nt] =gmsh.model.mesh.getNode(printNode)


    gmsh.finalize()
    macro.sparseSetUp()

    macro.setPropFunction(number=1,function=lambda t: t,tmin=0,tmax=1000)
    macro.setDt(dt=0.1)



    for i in range(10):
        macro.timeincr()
        macro.newton(refResidual=1e-9)

    sol=macro.getSolution(geo.vertexType(),nt,1)

    return sol    


import unittest

class TestPlasticBrick(unittest.TestCase):
    def __init__(self, methodName: str = "PlasticCantilever"):
        self.places=8
        super().__init__(methodName)

    def test_plastic(self):
        res = run(1)
        
        self.assertAlmostEqual(res[0], -6.657005328226617, places=self.places, msg="u_x is wrong")
        self.assertAlmostEqual(res[1], 94.5332303256969, places=self.places, msg="u_y is wrong")
        self.assertAlmostEqual(res[2], -0.0007299404182032974, places=self.places, msg="u_z is wrong")
    def test_plastic_order2(self):
        res = run(2)
        
        self.assertAlmostEqual(res[0], -7.13446527883948, places=self.places, msg="u_x is wrong")
        self.assertAlmostEqual(res[1], 101.35890651733864, places=self.places, msg="u_y is wrong")
        self.assertAlmostEqual(res[2], -0.00030262933371412025, places=self.places, msg="u_z is wrong")


if __name__ == '__main__':
    unittest.main()