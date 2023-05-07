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



box = gmsh.model.occ.addBox(0, 0, 0, L, b, b)
surfloop, surfaces = gmsh.model.occ.getSurfaceLoops(box)

print(surfaces)



gmsh.model.occ.synchronize()
# for i in surfaces[0]:

#     clt, ct = gmsh.model.occ.getCurveLoops(i)
#     for j in ct[0]:
#         print(j)
#         gmsh.model.mesh.setTransfiniteCurve(j, 4)
#     gmsh.model.mesh.setTransfiniteSurface(i)


# gmsh.model.mesh.setTransfiniteVolume(box)


order = 1

gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
#gmsh.option.setNumber("Mesh.RecombineAll", 2)
gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2)
gmsh.option.setNumber("Mesh.ElementOrder", order)
gmsh.model.mesh.generate(3)

gmsh.fltk().run()



pathname = os.path.dirname(sys.argv[0])
currPath = os.path.abspath(pathname)
fesys= HierAMuS.FEMPy(pathname, "firstvolume")
fesys.setSolutionState()
fesys.setSolver(3)


gm = fesys.getMeshCommands().getFromGMESH()
gm.addVertices(gmsh)
gm.addBrickFiniteElements(gmsh, order, tag=box, material=1)

mesh = fesys.getMeshCommands()
mesh.getElementFormulations().addEL300_3DSolid(1, meshiddisp=1, disporder=1, mode=1)
mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1, E=100, nu=0.3)
mesh.addMaterial(matNum=1, matFormNum=1, elemFormNum=1)


mesh.setDegreesOfFreedom()

geofe = fesys.getMeshCommands().getGeometryCommands()
bcfe = fesys.getMeshCommands().getBoundaryConditions()

bounFaces = gm.getFaceNumbers(gmsh, tag=surfaces[0][1], ftype=4, order=1)
bcfe.singleBC(geofe.faceType(), number=bounFaces, meshId=1, dofs=[1,1,1], shapeOrder=1)
loadFaces = gm.getFaceNumbers(gmsh, tag=surfaces[0][0], ftype=4, order=1)
bcfe.singleLoad(geofe.faceType(), number=loadFaces, meshId=1, load=[0,1,0], propnum=1, shapeorder=1)


gmsh.finalize()
fesys.getMacroCommands().sparseSetUp()

fesys.getMacroCommands().setPropFunction(number=1)
fesys.getMacroCommands().setDt(1)
fesys.getMacroCommands().timeincr()

fesys.getMacroCommands().setLogLevel(fesys.FullLog(), fesys.FullLog())
#fesys.getMacroCommands().newton(maxIteration=20)
fesys.getMacroCommands().assembleSolve()


fesys.getPlotCommands().toFile()


#fesys.getMacroCommands().printInfo()

#gmsh.fltk.run()