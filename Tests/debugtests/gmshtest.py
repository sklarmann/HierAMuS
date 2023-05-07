# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

# ------------------------------------------------------------------------------
#
#  Gmsh Python tutorial 11
#
#  Unstructured quadrangular meshes
#
# ------------------------------------------------------------------------------

import gmsh
import sys, os
import HierAMuS


gmsh.initialize()

gmsh.model.add("t11")

# We have seen in tutorials `t3.py' and `t6.py' that extruded and transfinite
# meshes can be "recombined" into quads, prisms or hexahedra. Unstructured
# meshes can be recombined in the same way. Let's define a simple geometry with
# an analytical mesh size field:

esize = 0.5

L=10
h=1

disporder = 2


p1 = gmsh.model.geo.addPoint(0,0,0,esize)
p2 = gmsh.model.geo.addPoint(L/2,0,0,esize)
p3 = gmsh.model.geo.addPoint(L/2,h,0,esize/4)
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
pl = gmsh.model.geo.addPlaneSurface([cl])
p2 = gmsh.model.geo.addPlaneSurface([cl2])

gmsh.model.geo.synchronize()

gmsh.option.setNumber("Mesh.Algorithm", 8)
gmsh.option().setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)

#field = gmsh.model.mesh.field
#field.add("MathEval", 1)
#field.setString(1, "F", "0.01*(1.0+30.*(y-x*x)*(y-x*x) + (1-x)*(1-x))")
#field.setAsBackgroundMesh(1)

# To generate quadrangles instead of triangles, we can simply add
#gmsh.model.mesh.setRecombine(2, pl)

# If we'd had several surfaces, we could have used the global option
# "Mesh.RecombineAll":
#
# gmsh.option.setNumber("Mesh.RecombineAll", 1)

# The default recombination algorithm is called "Blossom": it uses a minimum
# cost perfect matching algorithm to generate fully quadrilateral meshes from
# triangulations. More details about the algorithm can be found in the
# following paper: J.-F. Remacle, J. Lambrechts, B. Seny, E. Marchandise,
# A. Johnen and C. Geuzaine, "Blossom-Quad: a non-uniform quadrilateral mesh
# generator using a minimum cost perfect matching algorithm", International
# Journal for Numerical Methods in Engineering 89, pp. 1102-1119, 2012.

# For even better 2D (planar) quadrilateral meshes, you can try the
# experimental "Frontal-Delaunay for quads" meshing algorithm, which is a
# triangulation algorithm that enables to create right triangles almost
# everywhere: J.-F. Remacle, F. Henrotte, T. Carrier-Baudouin, E. Bechet,
# E. Marchandise, C. Geuzaine and T. Mouton. A frontal Delaunay quad mesh
# generator using the L^inf norm. International Journal for Numerical Methods
# in Engineering, 94, pp. 494-512, 2013. Uncomment the following line to try
# the Frontal-Delaunay algorithms for quads:
#
# gmsh.option.setNumber("Mesh.Algorithm", 8)

# The default recombination algorithm might leave some triangles in the mesh, if
# recombining all the triangles leads to badly shaped quads. In such cases, to
# generate full-quad meshes, you can either subdivide the resulting hybrid mesh
# (with `Mesh.SubdivisionAlgorithm' set to 1), or use the full-quad
# recombination algorithm, which will automatically perform a coarser mesh
# followed by recombination, smoothing and subdivision. Uncomment the following
# line to try the full-quad algorithm:
#
# gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2) # or 3

# You can also set the subdivision step alone, with
#
# gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)

gmsh.model.mesh.generate(2)


elementType = gmsh.model.mesh.getElementType("Quadrangle", 1)


pathname = os.path.dirname(sys.argv[0])

currPath = os.path.abspath(pathname)

fesys = HierAMuS.FEMPy(currPath, "gmshtest")

fesys.getMacroCommands().setLogLevel(fesys.FullLog(), fesys.FullLog())
fesys.setStaticSolutionState()

gm = fesys.getMeshCommands().getFromGMESH()

gm.geomFromGMSH(gmsh)


eltype = gmsh.model.mesh().getElementType("Quadrangle", 1)
gmsh.model.mesh.getElementsByType(elementType)

gm.addFiniteElementsFromTag(gmsh, eltype, pl, 1)
gm.addFiniteElementsFromTag(gmsh, eltype, p2, 1)

fesys.getMeshCommands().getMaterialFormulations().addMA3_2D_LinearElastic_Isotrop(1, E=100, nu=0.0, thickness=1, plainstrain=0)
fesys.getMeshCommands().getElementFormulations().addEL201_2DShell(1, meshiddisp=1, disporder=disporder, mode=3)
fesys.getMeshCommands().addMaterial(1,1,1)


fesys.getMeshCommands().setDegreesOfFreedom()

bedges = gm.getEdgeNumbers(gmsh, l4)
geo = fesys.getMeshCommands().getGeometryCommands()
fesys.getMeshCommands().getBoundaryConditions().singleBC(geo.edgeType(), bedges, 1, [1,1,1], disporder)

bedges = gm.getEdgeNumbers(gmsh, l6)
fesys.getMeshCommands().getBoundaryConditions().singleLoad(geo.edgeType(),bedges,1,[0,1,0],0,shapeorder=disporder)


fesys.getMacroCommands().sparseSetUp()

fesys.getMacroCommands().setPropFunction(0)
fesys.getMacroCommands().setDt(1)
fesys.getMacroCommands().timeincr()

fesys.getMacroCommands().assembleSolve()

fesys.getPlotCommands().toFile()


#fesys.getMacroCommands().printInfo()

nt,b,c = gmsh.model.mesh.getNodes(1,l2,includeBoundary=True)
print(nt)


#gmsh.fltk.run()
#fesys.getMeshCommands().getGeometryCommands().addLinearEdgeGeo(number, vertexList)







gmsh.finalize()