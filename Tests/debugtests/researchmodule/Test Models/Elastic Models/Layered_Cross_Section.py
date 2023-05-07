import os,sys
import FEMPyBind
import gmsh
#import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt

gmsh.initialize()
gmsh.model.add("test")

# FE
order=2
E=100
nu=0.0

# model
L=1
h1=0.4
h2=0.2
h3=0.4
b=1

nx=1
ny=1
nz=1

nz2=1

##Volume1 Points
gmsh.model.occ.addPoint(-L/2,-b/2,-h1-h2/2)         #1
gmsh.model.occ.addPoint(L/2,-b/2,-h1-h2/2)          #2
gmsh.model.occ.addPoint(L/2,b/2,-h1-h2/2)           #3
gmsh.model.occ.addPoint(-L/2,b/2,-h1-h2/2)          #4
gmsh.model.occ.addPoint(-L/2,-b/2,-h2/2)            #5
gmsh.model.occ.addPoint(L/2,-b/2,-h2/2)             #6
gmsh.model.occ.addPoint(L/2,b/2,-h2/2)              #7
gmsh.model.occ.addPoint(-L/2,b/2,-h2/2)             #8
##Volume2 top Points
gmsh.model.occ.addPoint(-L/2,-b/2,h2/2)             #9
gmsh.model.occ.addPoint(L/2,-b/2,h2/2)              #10
gmsh.model.occ.addPoint(L/2,b/2,h2/2)               #11
gmsh.model.occ.addPoint(-L/2,b/2,h2/2)              #12
##Volume3 top Points
gmsh.model.occ.addPoint(-L/2,-b/2,h3+h2/2)          #13
gmsh.model.occ.addPoint(L/2,-b/2,h3+h2/2)           #14
gmsh.model.occ.addPoint(L/2,b/2,h3+h2/2)            #15
gmsh.model.occ.addPoint(-L/2,b/2,h3+h2/2)           #16

#gmsh.model.occ.synchronize()                       ##check geometry##
#gmsh.fltk.run()                                    ##C.G. alligned with ref line and shear center##
                                                    ##check the material tangent at the end, should be diagonal## 
##Volume1 Faces
# Lines of lower face 1
gmsh.model.occ.addLine(1,2)
gmsh.model.occ.addLine(2,3)
gmsh.model.occ.addLine(3,4)
gmsh.model.occ.addLine(4,1)

lowerFaceLines = [1,4,3,2]

# left face lines 2
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
##Volume2 Faces
# left face2 lines 7
gmsh.model.occ.addLine(5,9)
gmsh.model.occ.addLine(9,12)
gmsh.model.occ.addLine(12,8)

leftFaceLines2 = [13,14,15,6]

# right face2 lines 8
gmsh.model.occ.addLine(6,10)
gmsh.model.occ.addLine(10,11)
gmsh.model.occ.addLine(11,7)

rightFaceLines2 = [16,17,18,9]
frontFaceLines2 = [11,16,19,13]
backFaceLines2 = [12,18,20,15]

# top face2 lines 9
gmsh.model.occ.addLine(9,10)
gmsh.model.occ.addLine(11,12)

topFaceLines2 = [20,14,19,17]

##Volume3 Faces
# left face3 lines 10
gmsh.model.occ.addLine(9,13)
gmsh.model.occ.addLine(13,16)
gmsh.model.occ.addLine(16,12)

leftFaceLines3 = [21,22,23,14]

# right face3 lines 11
gmsh.model.occ.addLine(10,14)
gmsh.model.occ.addLine(14,15)
gmsh.model.occ.addLine(15,11)

rightFaceLines3 = [24,25,26,17]
frontFaceLines3 = [19,24,27,21]
backFaceLines3 = [20,26,28,23]

# top face3 lines 12
gmsh.model.occ.addLine(13,14)
gmsh.model.occ.addLine(15,16)

topFaceLines3 = [28,22,27,25]
##End Volumes Faces
##Lines Lists
xlines = [1,3,11,12]
xlines2 = [11,12,19,20]
xlines3 = [19,20,27,28]

ylines = [2,9,4,6]
ylines2 = [9,17,6,14]
ylines3 = [17,14,25,22]

zlines = [5,7,8,10]
zlines2 = [13,15,16,18]
zlines3 = [21,23,24,26]
##End Lines Lists
##Volume 1
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

##Volume 2
# lower face2 1
c2=gmsh.model.occ.addCurveLoop(topFaceLines)
f2=gmsh.model.occ.addPlaneSurface([c2])

# left face2 2
c2=gmsh.model.occ.addCurveLoop(leftFaceLines2)
bounFace=gmsh.model.occ.addPlaneSurface([c2])

# right face2 3
c2=gmsh.model.occ.addCurveLoop(rightFaceLines2)
loadFace=gmsh.model.occ.addPlaneSurface([c2])

# top face2 4
c2=gmsh.model.occ.addCurveLoop(topFaceLines2)
gmsh.model.occ.addPlaneSurface([c2])

# front face2 5
c2=gmsh.model.occ.addCurveLoop(frontFaceLines2)
gmsh.model.occ.addPlaneSurface([c2])
# top face2 6
c2=gmsh.model.occ.addCurveLoop(backFaceLines2)
gmsh.model.occ.addPlaneSurface([c2])

gmsh.model.occ.addSurfaceLoop([7,8,9,10,11,12])
gmsh.model.occ.addVolume([2])

##Volume 3
# lower face2 1
c3=gmsh.model.occ.addCurveLoop(topFaceLines2)
f3=gmsh.model.occ.addPlaneSurface([c3])

# left face2 2
c3=gmsh.model.occ.addCurveLoop(leftFaceLines3)
bounFace=gmsh.model.occ.addPlaneSurface([c3])

# right face2 3
c3=gmsh.model.occ.addCurveLoop(rightFaceLines3)
loadFace=gmsh.model.occ.addPlaneSurface([c3])

# top face2 4
c3=gmsh.model.occ.addCurveLoop(topFaceLines3)
gmsh.model.occ.addPlaneSurface([c3])

# front face2 5
c3=gmsh.model.occ.addCurveLoop(frontFaceLines3)
gmsh.model.occ.addPlaneSurface([c3])
# top face2 6
c3=gmsh.model.occ.addCurveLoop(backFaceLines3)
gmsh.model.occ.addPlaneSurface([c3])

gmsh.model.occ.addSurfaceLoop([13,14,15,16,17,18])
gmsh.model.occ.addVolume([3])

##End volumes definition

gmsh.model.occ.synchronize()


# setting the lines to be transfinite Vol1
for i in xlines:
    gmsh.model.mesh.setTransfiniteCurve(i,nx+1)
for i in ylines:
    gmsh.model.mesh.setTransfiniteCurve(i,ny+1)
for i in zlines:
    gmsh.model.mesh.setTransfiniteCurve(i,nz+1)

# setting the lines to be transfinite Vol2       ##vol 2 refinement should be independent on the other volumes in z dir#
for i in xlines2:
    gmsh.model.mesh.setTransfiniteCurve(i,nx+1)
for i in ylines2:
    gmsh.model.mesh.setTransfiniteCurve(i,ny+1)
for i in zlines2:
    gmsh.model.mesh.setTransfiniteCurve(i,nz2+1)

# setting the lines to be transfinite Vol3
for i in xlines3:
    gmsh.model.mesh.setTransfiniteCurve(i,nx+1)
for i in ylines3:
    gmsh.model.mesh.setTransfiniteCurve(i,ny+1)
for i in zlines3:
    gmsh.model.mesh.setTransfiniteCurve(i,nz+1)
    
# setting the faces to be transfinite
for i in [1,2,3,4,5,6]:
    gmsh.model.mesh.setTransfiniteSurface(i)
# setting the faces to be transfinite
for i in [7,8,9,10,11,12]:
    gmsh.model.mesh.setTransfiniteSurface(i)
# setting the faces to be transfinite
for i in [13,14,15,16,17,18]:
    gmsh.model.mesh.setTransfiniteSurface(i)
    
# setting the volume to be transfinite
gmsh.model.mesh.setTransfiniteVolume(1)
gmsh.model.mesh.setTransfiniteVolume(2)
gmsh.model.mesh.setTransfiniteVolume(3)


gmsh.option.setNumber("Mesh.RecombineAll", 1) # This option sets gmsh to recombine tetra to bricks

gmsh.model.mesh.generate(3)
gmsh.fltk.run()

currPath = os.path.dirname(os.path.realpath(__file__))

Rve = FEMPyBind.FEMPy(currPath,"RVE")
Rve.setStaticHomogenizationSolutionState()
Rve.setSolver(6)
Rve.getMacroCommands().setLogLevel(Rve.FullLog(),Rve.BasicLog())

mesh = Rve.getMeshCommands()
macro = Rve.getMacroCommands()
gm = mesh.getFromGMESH()
geo = mesh.getGeometryCommands()


gm.addGeomFromGmsh(gmsh)
geo.checkGeometry()

gm.addVolumeElements(gmsh,1,1)
gm.addVolumeElements(gmsh,2,2)
gm.addVolumeElements(gmsh,3,1)


mesh.getElementFormulations().addEL300_3DSolid(num=1,meshiddisp=1,disporder=order,mode=1)
mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1,E=E,nu=nu)
mesh.addMaterial(matNum=1,matFormNum=1,elemFormNum=1)

mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(2,E=E/10000,nu=nu)
mesh.addMaterial(matNum=2,matFormNum=2,elemFormNum=1)


mesh.setDegreesOfFreedom()

macro.getHomogenizationCommands().setHomogenizationBeam(meshIdDisp=1,dispOrder=order,bctype=0)

macro.sparseSetUp()

macro.getHomogenizationCommands().computeAMatrix()
macro.getHomogenizationCommands().setStrains([0.0,0,0,0,0.1,0])
eps=0
d_eps=0.1

macro.newton(refResidual=1e-11)

Rve.getPlotCommands().toFile()

macro.getHomogenizationCommands().homogenize()
macro.printInfo()