# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import HierAMuS
import gmsh
import os, sys

gmsh.initialize()
gmsh.option.setNumber('General.Terminal', 0)

gmsh.model.add("test")

# FE
order=2
E=100
nu=0.0

# model
L=10
h=1
b=1

nx=100
ny=10
nz=10




gmsh.model.occ.addPoint(-L/2,-b/2,-h/2)
gmsh.model.occ.addPoint(L/2,-b/2,-h/2)
gmsh.model.occ.addPoint(L/2,b/2,-h/2)
gmsh.model.occ.addPoint(-L/2,b/2,-h/2)
gmsh.model.occ.addPoint(-L/2,-b/2,h/2)
gmsh.model.occ.addPoint(L/2,-b/2,h/2)
gmsh.model.occ.addPoint(L/2,b/2,h/2)
gmsh.model.occ.addPoint(-L/2,b/2,h/2)

pf = gmsh.model.occ.addPoint(L/2,0,0)

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


# settings the lines to be transfinite
for i in xlines:
    gmsh.model.mesh.setTransfiniteCurve(i,nx+1)
for i in ylines:
    gmsh.model.mesh.setTransfiniteCurve(i,ny+1)
for i in zlines:
    gmsh.model.mesh.setTransfiniteCurve(i,nz+1)
    
# setting the faces to be transfinite
for i in [1,2,3,4,5,6]:
    gmsh.model.mesh.setTransfiniteSurface(i)
    
# setting the volume to be transfinite
gmsh.model.mesh.setTransfiniteVolume(1)


gmsh.option.setNumber("Mesh.RecombineAll", 1) # This option sets gmsh to recombine tetra to bricks

gmsh.model.mesh.generate(3)
#gmsh.fltk.run()


currPath = os.path.dirname(os.path.realpath(__file__))
print(currPath)
fesys = HierAMuS.FEMPy(currPath,"faceConstraint")
fesys.setStaticSolutionState()
fesys.setSolver(4)

macro = fesys.getMacroCommands()
mesh = fesys.getMeshCommands()
geo = mesh.getGeometryCommands()
gm = mesh.getFromGMESH()

gm.addGeomFromGmsh(gmsh)
geo.checkGeometry()

gm.addVolumeElements(gmsh,1,1)

gm.addFaceConstraint(gmsh,faceTags=loadFace,vertexTag=pf,material=2)


mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1,E,nu)
mesh.getElementFormulations().addEL300_3DSolid(1,1,order,1)
mesh.addMaterial(1,1,1)


mesh.getElementFormulations().addEL207_FaceConstraint(2,1,1,2,3,4,order,1)
mesh.addMaterial(2,1,2)

mesh.setDegreesOfFreedom()

facelist = gm.getFaceNumbers(gmsh,bounFace,4,1)
mesh.getBoundaryConditions().FaceBC(facelist,1,[1,1,1],order)

#mesh.getBoundaryConditions().BC(geo.vertexType(),pf,2,[1,1,1],1)
#mesh.getBoundaryConditions().BC(geo.vertexType(),pf,4,[1,1,1],1)

mesh.getBoundaryConditions().load(geo.vertexType(),pf,2,[1,0,0],0)

macro.sparseSetUp()

macro.setLogLevel(fesys.FullLog(),fesys.BasicLog())
macro.setPropFunction(0)
macro.setDt(1)

for i in range(0,10):
    macro.timeincr()
    macro.newton()

#macro.computeEigenValues(6,20)

fesys.getPlotCommands().toFile()

macro.printInfo()

