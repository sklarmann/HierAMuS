import os,sys
import FEMPyBind
import gmsh
import math
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
L=3
h=3
b=1.5
r=0.65
a=math.pi/4

xr=r*math.cos(a)
yr=r*math.sin(a)

nx=5
ny=5
nz=5

##Volume1 Points
gmsh.model.occ.addPoint(-L/2,-b/2,-h/2)             #1
gmsh.model.occ.addPoint(L/2,-b/2,-h/2)              #2
gmsh.model.occ.addPoint(-L/2,-b/2,h/2)              #3
gmsh.model.occ.addPoint(L/2,-b/2,h/2)               #4
gmsh.model.occ.addPoint(0,-b/2,0)                   #5
gmsh.model.occ.addPoint(-xr,-b/2,-yr)               #6
gmsh.model.occ.addPoint(xr,-b/2,-yr)                #7
gmsh.model.occ.addPoint(-xr,-b/2,yr)                #8
gmsh.model.occ.addPoint(xr,-b/2,yr)                 #9
gmsh.model.occ.addPoint(0,-b/2,-r)                  #10
gmsh.model.occ.addPoint(0,-b/2,r)                   #11
gmsh.model.occ.addPoint(-r,-b/2,0)                  #12
gmsh.model.occ.addPoint(r,-b/2,0)                   #13
gmsh.model.occ.addPoint(0,-b/2,-h/2)                #14
gmsh.model.occ.addPoint(0,-b/2,h/2)                 #15
gmsh.model.occ.addPoint(-L/2,-b/2,0)                #16
gmsh.model.occ.addPoint(L/2,-b/2,0)                 #17
#gmsh.model.occ.synchronize()                       ##check geometry##
#gmsh.fltk.run()  

##Lines
gmsh.model.occ.addLine(1,14)                       #L1
gmsh.model.occ.addLine(2,14)                       #L2
gmsh.model.occ.addLine(1,16)                       #L3
gmsh.model.occ.addLine(16,3)                       #L4
gmsh.model.occ.addLine(3,15)                       #L5
gmsh.model.occ.addCircleArc(8,5,11)                #L6
gmsh.model.occ.addCircleArc(11,5,9)                #L7
gmsh.model.occ.addCircleArc(9,5,13)                #L8
gmsh.model.occ.addCircleArc(13,5,7)                #L9
gmsh.model.occ.addCircleArc(7,5,10)                #L10
gmsh.model.occ.addCircleArc(10,5,6)                #L11
gmsh.model.occ.addCircleArc(6,5,12)                #L12
gmsh.model.occ.addCircleArc(12,5,8)                #L13
gmsh.model.occ.addLine(1,6)                        #L14
gmsh.model.occ.addLine(2,7)                        #L15
gmsh.model.occ.addLine(3,8)                        #L16
gmsh.model.occ.addLine(4,9)                        #L17
gmsh.model.occ.addLine(12,16)                      #L18
gmsh.model.occ.addLine(13,17)                      #L19
gmsh.model.occ.addLine(11,15)                      #L20
gmsh.model.occ.addLine(10,14)                      #L21
gmsh.model.occ.addLine(2,17)                       #L22
gmsh.model.occ.addLine(17,4)                       #L23
gmsh.model.occ.addLine(15,4)                       #24
gmsh.model.occ.synchronize()                       ##check geometry##
#gmsh.fltk.run()  

# front face lines 
frontFaceLines1 = [1,14,11,21]
frontFaceLines2 = [21,10,15,2]
frontFaceLines3 = [14,3,18,12]
frontFaceLines4 = [9,19,22,15]
frontFaceLines5 = [18,4,16,13]
frontFaceLines6 = [8,17,23,19]
frontFaceLines7 = [16,5,20,6]
frontFaceLines8 = [20,24,17,7]

gmsh.model.occ.synchronize()  
#gmsh.fltk.run()  

##Volume 1

c1=gmsh.model.occ.addCurveLoop(frontFaceLines1)
gmsh.model.occ.addPlaneSurface([c1])

c1=gmsh.model.occ.addCurveLoop(frontFaceLines2)
gmsh.model.occ.addPlaneSurface([c1])

c1=gmsh.model.occ.addCurveLoop(frontFaceLines3)
gmsh.model.occ.addPlaneSurface([c1])

c1=gmsh.model.occ.addCurveLoop(frontFaceLines4)
gmsh.model.occ.addPlaneSurface([c1])

c1=gmsh.model.occ.addCurveLoop(frontFaceLines5)
gmsh.model.occ.addPlaneSurface([c1])

c1=gmsh.model.occ.addCurveLoop(frontFaceLines6)
gmsh.model.occ.addPlaneSurface([c1])

c1=gmsh.model.occ.addCurveLoop(frontFaceLines7)
gmsh.model.occ.addPlaneSurface([c1])

c1=gmsh.model.occ.addCurveLoop(frontFaceLines8)
gmsh.model.occ.addPlaneSurface([c1])
gmsh.model.occ.synchronize()

# setting the faces to be transfinite before exruding
for i in range(24):
    print(i+1)
    gmsh.model.mesh.setTransfiniteCurve(i+1,4)
#gmsh.fltk.run()
for i in [1,2,3,4,5,6,7,8]:
    gmsh.model.mesh.setTransfiniteSurface(i)
gmsh.model.occ.synchronize() 
gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
gmsh.model.occ.synchronize() 
gmsh.option.setNumber("Mesh.RecombineAll", 1) # This option sets gmsh to recombine tetra to bricks
gmsh.model.occ.synchronize()  
#gmsh.fltk.run()

# extruding faces 
gmsh.model.occ.extrude([(2,1),(2,2),(2,3),(2,4),(2,5),(2,6),(2,7),(2,8)],0,3,0,[4],recombine=True)
gmsh.model.occ.synchronize()  
#gmsh.fltk.run()
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()  
gmsh.fltk.run()

gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.RecombineAll", 1) # This option sets gmsh to recombine tetra to bricks
gmsh.model.occ.synchronize()  
#Generate 3D Mesh
gmsh.model.mesh.generate(3)
gmsh.model.occ.synchronize()  
gmsh.fltk.run()
##############################################################################################
# Copying volumes
# gmsh.model.occ.copy([(2,1),(2,2),(2,3),(2,4),(2,5),(2,6),(2,7),(2,8)])
# gmsh.model.occ.translate([(3,1),(3,2),(3,3),(3,4),(3,5),(3,6),(3,7),(3,8)],0,2,0)
# gmsh.model.occ.synchronize()  
#####
# for i in [41,42,43,44,45,46,47,48]:
#     gmsh.model.mesh.setTransfiniteSurface(i)
# gmsh.model.occ.synchronize() 
# gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
# gmsh.model.occ.synchronize() 
# gmsh.option.setNumber("Mesh.RecombineAll", 1) # This option sets gmsh to recombine tetra to bricks
# gmsh.model.occ.synchronize()  
# gmsh.fltk.run()
# gmsh.model.occ.extrude([(3,9),(3,10),(3,11),(3,12),(3,13),(3,14),(3,15),(3,16)],0,.2,0,[1],recombine=True)
# gmsh.model.occ.synchronize()  
# gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
# gmsh.model.occ.synchronize()
# gmsh.option.setNumber("Mesh.RecombineAll", 1) # This option sets gmsh to recombine tetra to bricks

##### 
#gmsh.model.occ.synchronize()  
# # setting the faces to be transfinite before exruding
# for i in [85,90,95,100,105,110,115,119]:
#     gmsh.model.mesh.setTransfiniteSurface(i)
# gmsh.model.occ.synchronize() 
# gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
# gmsh.model.occ.synchronize() 
# gmsh.option.setNumber("Mesh.RecombineAll", 1) # This option sets gmsh to recombine tetra to bricks
# gmsh.model.occ.synchronize()  
# gmsh.model.mesh.generate(3)
# gmsh.model.occ.synchronize()  
#gmsh.fltk.run()

# setting the faces to be transfinite after exruding
# for i in [9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]:
#     gmsh.model.mesh.setTransfiniteSurface(i)
# gmsh.model.occ.synchronize() 

# gmsh.model.mesh.generate(3)
# gmsh.model.occ.synchronize()  
##############################################################################################

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
Vols=[(3,1),(3,2),(3,3),(3,4),(3,5),(3,6),(3,7),(3,8)]

for x_i in range(1, len(Vols)+1, 1):
    gm.addVolumeElements(gmsh,x_i,1)


################################################## uncomment later

mesh.getElementFormulations().addEL300_3DSolid(num=1,meshiddisp=1,disporder=order,mode=1)
mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1,E=E,nu=nu)
mesh.addMaterial(matNum=1,matFormNum=1,elemFormNum=1)

mesh.setDegreesOfFreedom()

macro.getHomogenizationCommands().setHomogenizationSolid3D(meshIdDisp=1,dispOrder=order,bctype=0)

macro.sparseSetUp()

macro.getHomogenizationCommands().computeAMatrix()
macro.getHomogenizationCommands().setStrains([0.1,0,0,0,0,0])
eps=0
d_eps=0.1

macro.newton(refResidual=1e-11)

Rve.getPlotCommands().toFile()

macro.getHomogenizationCommands().homogenize()
macro.printInfo()