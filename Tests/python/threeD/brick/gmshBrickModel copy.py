# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import gmsh
import HierAMuS
import os, sys

gmsh.initialize()

gmsh.model.add("t11")

gmsh.model.setFileName(fileName="/home/simon/github/cppfemtemp/Tests/python/3D/brick/t11")

L=10
b=1
h=1
fac=4
nx = fac*10
ny = fac*2
nz = fac*2

points = []



points.append(gmsh.model.occ.addPoint(0, 0, 0))
points.append(gmsh.model.occ.addPoint(L, 0, 0))
points.append(gmsh.model.occ.addPoint(L, b, 0))
points.append(gmsh.model.occ.addPoint(0, b, 0))
points.append(gmsh.model.occ.addPoint(0, 0, h))
points.append(gmsh.model.occ.addPoint(L, 0, h))
points.append(gmsh.model.occ.addPoint(L, b, h))
points.append(gmsh.model.occ.addPoint(0, b, h))

blines = []
tlines = []
clines = []

blines.append(gmsh.model.occ.addLine(points[0], points[1]))
blines.append(gmsh.model.occ.addLine(points[1], points[2]))
blines.append(gmsh.model.occ.addLine(points[2], points[3]))
blines.append(gmsh.model.occ.addLine(points[3], points[0]))

tlines.append(gmsh.model.occ.addLine(points[4], points[5]))
tlines.append(gmsh.model.occ.addLine(points[5], points[6]))
tlines.append(gmsh.model.occ.addLine(points[6], points[7]))
tlines.append(gmsh.model.occ.addLine(points[7], points[4]))

clines.append(gmsh.model.occ.addLine(points[0], points[4]))
clines.append(gmsh.model.occ.addLine(points[1], points[5]))
clines.append(gmsh.model.occ.addLine(points[2], points[6]))
clines.append(gmsh.model.occ.addLine(points[3], points[7]))



curveLoops = []
curveLoops.append(gmsh.model().occ.addCurveLoop(blines)) # bottom
curveLoops.append(gmsh.model().occ.addCurveLoop(tlines)) # top
curveLoops.append(gmsh.model().occ.addCurveLoop([blines[0],clines[1],-clines[0],-tlines[0]]))
curveLoops.append(gmsh.model().occ.addCurveLoop([blines[1],clines[2],-clines[1],-tlines[1]])) #front
curveLoops.append(gmsh.model().occ.addCurveLoop([blines[2],clines[3],-clines[2],-tlines[2]]))
curveLoops.append(gmsh.model().occ.addCurveLoop([blines[3],clines[0],-clines[3],-tlines[3]])) # back
faces = []
for i in curveLoops:
    faces.append(gmsh.model.occ.addPlaneSurface([i]))

fl=gmsh.model.occ.addSurfaceLoop(faces)
vol=gmsh.model.occ.addVolume([fl])

gmsh.model.occ.synchronize()

# lines in length direction
llines = [blines[0],blines[2],tlines[0],tlines[2]]
for l in llines:
    gmsh.model.mesh.setTransfiniteCurve(l, nx+1)

# lines in width direction
llines = [blines[1],blines[3],tlines[1],tlines[3]]
for l in llines:
    gmsh.model.mesh.setTransfiniteCurve(l, ny+1)

# lines in height direction
llines = clines
for l in llines:
    gmsh.model.mesh.setTransfiniteCurve(l, nz+1)



for f in faces:
    gmsh.model.mesh.setTransfiniteSurface(f)

gmsh.model.mesh.setTransfiniteVolume(vol)

order = 1

gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
gmsh.option.setNumber("Mesh.RecombineAll", 2)
#gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2)
gmsh.option.setNumber("Mesh.ElementOrder", order)
gmsh.model.mesh.generate(3)




pathname = os.path.dirname(sys.argv[0])
currPath = os.path.abspath(pathname)
fesys= HierAMuS.FEMPy(pathname, "firstvolume")
fesys.setSolutionState()
fesys.setSolver(3)


gm = fesys.getMeshCommands().getFromGMESH()
gm.addVertices(gmsh)
gm.addBrickFiniteElements(gmsh, order, tag=vol, material=1)

mesh = fesys.getMeshCommands()
mesh.getElementFormulations().addEL300_3DSolid(1, meshiddisp=1, disporder=1, mode=2)
mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1, E=100, nu=0.3)
mesh.addMaterial(matNum=1, matFormNum=1, elemFormNum=1)


mesh.setDegreesOfFreedom()

geofe = fesys.getMeshCommands().getGeometryCommands()
bcfe = fesys.getMeshCommands().getBoundaryConditions()

bounFaces = gm.getFaceNumbers(gmsh, tag=faces[5], ftype=4, order=1)
bcfe.singleBC(geofe.faceType(), number=bounFaces, meshId=1, dofs=[1,1,1], shapeOrder=1)
loadFaces = gm.getFaceNumbers(gmsh, tag=faces[3], ftype=4, order=1)
bcfe.singleLoad(geofe.faceType(), number=loadFaces, meshId=1, load=[0,1,0], propnum=1, shapeorder=1)


gmsh.finalize()
fesys.getMacroCommands().sparseSetUp()

fesys.getMacroCommands().setPropFunction(number=1)
fesys.getMacroCommands().setDt(1)
fesys.getMacroCommands().timeincr()

fesys.getMacroCommands().setLogLevel(fesys.FullLog(), fesys.FullLog())
fesys.getMacroCommands().newton(maxIteration=20)

fesys.getPlotCommands().toFile()


#fesys.getMacroCommands().printInfo()

#gmsh.fltk.run()