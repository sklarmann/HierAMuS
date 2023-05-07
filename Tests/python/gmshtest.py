# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import gmsh
import sys, os
import HierAMuS
import numpy as np


gmsh.initialize()

gmsh.model.add("t11")

# We have seen in tutorials `t3.py' and `t6.py' that extruded and transfinite
# meshes can be "recombined" into quads, prisms or hexahedra. Unstructured
# meshes can be recombined in the same way. Let's define a simple geometry with
# an analytical mesh size field:

esize = 0.1

L=2
h=1

disporder = 4
order=1


p1 = gmsh.model.geo.addPoint(0,0,0,esize)
p2 = gmsh.model.geo.addPoint(L/2,0,0,esize)
p3 = gmsh.model.geo.addPoint(L/2,h,0,esize)
p4 = gmsh.model.geo.addPoint(0,h,0,esize)



p5 = gmsh.model.geo.addPoint(L,0,0,esize)
p6 = gmsh.model.geo.addPoint(L,h,0,esize)


l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)

l5 = gmsh.model.geo.addLine(p2, p5)
l6 = gmsh.model.geo.addLine(p5, p6)
l7 = gmsh.model.geo.addLine(p6, p3)

cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
cl2 = gmsh.model.geo.addCurveLoop([l5, l6, l7, -l2])
p1 = gmsh.model.geo.addPlaneSurface([cl])
p2 = gmsh.model.geo.addPlaneSurface([cl2])

gmsh.model.geo.synchronize()

gmsh.option.setNumber("Mesh.Algorithm", 11)
gmsh.option().setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
gmsh.option.setNumber("Mesh.ElementOrder", order)

gmsh.model.mesh.generate(2)


pathname = os.path.dirname(sys.argv[0])
currPath = os.path.abspath(pathname)
fesys= HierAMuS.FEMPy(pathname, "firstvolume")
fesys.setSolutionState()

gm = fesys.getMeshCommands().getFromGMESH()
gm.addVertices(gmsh)
gm.addQuadrilateralFiniteElements(gmsh, order, p1, 1)
gm.addQuadrilateralFiniteElements(gmsh, order, p2, 1)

fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(number=1, E=100, nu=0.3, thickness=1, plainstrain=0)
fesys.getMeshCommands().getElementFormulations().addEL201_2DShell(num=1, meshiddisp=1, disporder=disporder, mode=1)
fesys.getMeshCommands().addMaterial(matNum=1, matFormNum=1, elemFormNum=1)

fesys.getMeshCommands().setDegreesOfFreedom()

enums = gm.getEdgeNumbers(gmsh, tag=l4, order=order)
fesys.getMeshCommands().getBoundaryConditions().singleBC(fesys.getMeshCommands().getGeometryCommands().edgeType(), number=enums, meshId=1, dofs=[1,1,1], shapeOrder=disporder)
enums = gm.getEdgeNumbers(gmsh, tag=l6, order=order)
fesys.getMeshCommands().getBoundaryConditions().singleLoad(fesys.getMeshCommands().getGeometryCommands().edgeType(), number=enums, meshId=1, load=[0,1,0], propnum=0,shapeorder=disporder)
fesys.getMacroCommands().sparseSetUp()


fesys.getMacroCommands().setPropFunction(0)
fesys.getMacroCommands().setDt(1)
fesys.getMacroCommands().timeincr()

fesys.getMacroCommands().assembleSolve()

fesys.getPlotCommands().toFile()
gmsh.model.mesh.getAllEdges()


# gmsh.model.mesh.createEdges()
# eltype = gmsh.model.mesh.getElementType("Quadrangle", order)
# facetags, facenodes = gmsh.model.mesh.getElementsByType(eltype,tag=p2)
# numFaceNodes = (order+1)**2
# print(facetags)
# for i in range(len(facetags)):
#     fnum = facetags[i]
#     print(fnum)
#     fn = facenodes[i*numFaceNodes:(i+1)*numFaceNodes]
#     el = fn[0:2]
#     el = np.append(el,fn[1:3])
#     el = np.append(el,fn[2:4])
#     el = np.append(el,[fn[3],fn[0]])
#     et,eo = gmsh.model.mesh.getEdges(el)
#     for i in range(len(et)):
#         enum = et[i]
#         everts = el[2*i:2*i+2]
#         fesys.getMeshCommands().getGeometryCommands().addLinearEdgeGeo(enum, everts)
#         print(everts)
#     fesys.getMeshCommands().getGeometryCommands().addQuadrilateralFace(fnum, fn, et)




fesys.getMacroCommands().setLogLevel(fesys.FullLog(), fesys.FullLog())
fesys.getMacroCommands().printInfo()
gmsh.fltk().run()
