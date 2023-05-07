# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import os,sys
import HierAMuS
import gmsh

currPath = os.path.dirname(os.path.realpath(__file__))


gmsh.initialize()
gmsh.model.add("test")

# FE
order=1
E=100.0
nu=0.0
#alpha = 22.8*10**-6
alpha=1.0
#alpha=0
#c = 452
c = 0
rho0 = 0
T0 = 0.0
kappa = 1.0

lam = E*nu/((1+nu)*(1-2*nu))
G=E/(2.0*(1.0+nu))

# model
L=20
h=1
b=1
d=0.1

nx=100
ny=5
nz=5

gmsh.model.occ.addBox(0,-b/2,-h/2,L-d,b,h)
gmsh.model.occ.addBox(L-d,-b/2,-h/2,d,b,h)

gmsh.model.occ.removeAllDuplicates()


xlines1 = [9,10,11,12]
xlines2 = [17,18,19,20]
ylines = [2,4,6,8,14,16]
zlines = [1,3,5,7,13,15]

xfaces = [1,2,7]
bounface = 1
loadface = 7
yfaces1 = [5,6]
yfaces2 = [10,11]
zfaces1 = [3,4]
zfaces2 = [8,9]

gmsh.model.occ.synchronize()

for i in xlines1:
    gmsh.model.mesh.setTransfiniteCurve(i,nx)
for i in xlines2:
    gmsh.model.mesh.setTransfiniteCurve(i,2)

for i in ylines:
    gmsh.model.mesh.setTransfiniteCurve(i,ny)
for i in zlines:
    gmsh.model.mesh.setTransfiniteCurve(i,nz)
    
for i in xfaces:
    gmsh.model.mesh.setTransfiniteSurface(i)
for i in yfaces1:
    gmsh.model.mesh.setTransfiniteSurface(i)
for i in yfaces2:
    gmsh.model.mesh.setTransfiniteSurface(i)
for i in zfaces1:
    gmsh.model.mesh.setTransfiniteSurface(i)
for i in zfaces2:
    gmsh.model.mesh.setTransfiniteSurface(i)

gmsh.model.mesh.setTransfiniteVolume(1)
gmsh.model.mesh.setTransfiniteVolume(2)


gmsh.option.setNumber("Mesh.RecombineAll", 1) # This option sets gmsh to recombine tetra to bricks

gmsh.model.mesh.generate(3)

gmsh.fltk.run()

fesys = HierAMuS.FEMPy(currPath,"BeamModel3D2")
fesys.setStaticSolutionState()
fesys.setSolver(6)
fesys.getMacroCommands().setLogLevel(fesys.FullLog(),fesys.BasicLog())

mesh = fesys.getMeshCommands()
macro = fesys.getMacroCommands()
gm = mesh.getFromGMESH()
geo = mesh.getGeometryCommands()


gm.addGeomFromGmsh(gmsh)
geo.checkGeometry()



gm.addVolumeElements(gmsh,1,1)
gm.addVolumeElements(gmsh,2,1)
mesh.getElementFormulations().addEL303_ThermoMechanikSolid3D(num=1,meshiddisp=1,meshidtemperatur=2,shapeorder=order,mu=G,lamb=lam,alpha=alpha,c=c,rho0=rho0,T0=T0,kappa=kappa,mode=1)
mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(number=1,E=E,nu=nu)
mesh.addMaterial(matNum=1,matFormNum=1,elemFormNum=1)

mesh.setDegreesOfFreedom()

fnums = gm.getFaceNumbers(gmsh,bounface,4,1)
mesh.getBoundaryConditions().BCFace(fnums,1,[1,1,1],order)
mesh.getBoundaryConditions().BCFace(fnums,2,[1,1,1],order)

fnums = gm.getFaceNumbers(gmsh,loadface,4,1)
mesh.getBoundaryConditions().BCFace(fnums,1,[1,0,0],order)
mesh.getBoundaryConditions().BCFace(fnums,2,[1,1,1],order)
#mesh.getBoundaryConditions().PSFaceH1(fnums,2,[1,0,0],1,add=False)
#fnums = gm.getFaceNumbers(gmsh,zfaces2[0],4,1)
#mesh.getBoundaryConditions().BCFace(fnums,2,[1,1,1],order)
#mesh.getBoundaryConditions().PSFaceH1(fnums,2,[1,0,0],1,add=False)
#fnums = gm.getFaceNumbers(gmsh,zfaces2[1],4,1)
#mesh.getBoundaryConditions().BCFace(fnums,2,[1,1,1],order)
#mesh.getBoundaryConditions().PSFaceH1(fnums,2,[-1,0,0],1,add=False)

nt,coord,paramCoord=gmsh.model.mesh.getNodes(2,loadface,includeBoundary=True)

for i in range(len(nt)):
    vnum = nt[i]
    ycoor = coord[3*i+1]
    lval = 2*ycoor
    mesh.getBoundaryConditions().PSVertex(number=vnum,meshId=2,load=[-lval,0,0],propnum=1,add=False)
    

macro.sparseSetUp()



print(nt)
print(coord)

macro.setPropFunction(1)
macro.setDt(1)

macro.timeincr()

macro.newton(refResidual=1e-9)

fesys.getPlotCommands().toFile()

macro.printInfo()

E=100
I=1

F=1
print("ref",F*L**3/E/I/3+F*L/50)

print(macro.getSolution(geo.vertexType(),2,3))
