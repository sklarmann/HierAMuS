import os,sys
import HierAMuS
import gmsh

gmsh.initialize()
gmsh.model.add("test")

# FE
order=2
E=100
nu=0.0

# model
fac=100
L=fac
b=fac
h=1

nx=8
ny=8
nz=4




gmsh.model.occ.addPoint(-L/2,-b/2,-h/2)
gmsh.model.occ.addPoint(L/2,-b/2,-h/2)
gmsh.model.occ.addPoint(L/2,b/2,-h/2)
gmsh.model.occ.addPoint(-L/2,b/2,-h/2)
gmsh.model.occ.addPoint(-L/2,-b/2,h/2)
gmsh.model.occ.addPoint(L/2,-b/2,h/2)
gmsh.model.occ.addPoint(L/2,b/2,h/2)
gmsh.model.occ.addPoint(-L/2,b/2,h/2)

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

Rve = HierAMuS.FEMPy(currPath,"RVE")
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
mesh.getElementFormulations().addEL300_3DSolid(num=1,meshiddisp=1,disporder=order,mode=1)
mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1,E=E,nu=nu)
mesh.addMaterial(matNum=1,matFormNum=1,elemFormNum=1)

mesh.setDegreesOfFreedom()

macro.getHomogenizationCommands().setHomogenizationShell(meshIdDisp=1,dispOrder=order,bctype=0)

macro.sparseSetUp()

macro.getHomogenizationCommands().computeAMatrix()
macro.getHomogenizationCommands().setStrains([0,0,0,0,0.1,0,0,0])

macro.newton(refResidual=1e-11)

Rve.getPlotCommands().toFile()

macro.getHomogenizationCommands().homogenize()
macro.printInfo()