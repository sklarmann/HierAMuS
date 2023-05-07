# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from sympy import true
import HierAMuS
from HierAMuS.HierAMuSPyWrapper import HierAMuSPyFEM

# Hier ggf. zuvor "python -m pip install matplotlib" in cmd ausfuehren
import matplotlib.pyplot as plt
import matplotlib

import os, sys
import time

#time.sleep(20)

print('sys.argv[0] =', sys.argv[0])             
pathname = os.path.dirname(sys.argv[0])        
print('path =', pathname)
print('full path =', os.path.abspath(pathname)) 

currPath = os.path.abspath(pathname)


# length of the beam
L   = 10
# transition element length
#ell = 0.001
# layer heights
h=1  


nx = 2
ny = 2
sy = 0

fac=1

nx*=fac
ny*=fac


# load
meshid=1
order=2
load=[0,1,0]

# Anlegen des FE Systems mit Pfad und Dateiname, in die Daten geschrieben.
fesys = HierAMuS.FEMPy(currPath, 'test.log')
# Loesungstyp festlegen
fesys.setStaticSolutionState()


geoB = fesys.getGeomBuilder()

v1 = geoB.getVertex()
v2 = geoB.getVertex()
v3 = geoB.getVertex()
v4 = geoB.getVertex()
v1.setCoordinates(0,0,0)
v2.setCoordinates(L,0,0)
v3.setCoordinates(L,h,0)
v4.setCoordinates(0,h,0)

e1 = geoB.getEdge()
e2 = geoB.getEdge()
e3 = geoB.getEdge()
e4 = geoB.getEdge()
e1.setStartEnd(v1,v2)
e2.setStartEnd(v2,v3)
e3.setStartEnd(v3,v4)
e4.setStartEnd(v4,v1)


q1 = geoB.getQuad()
q1.setVertsEdges([v1,v2,v3,v4],[e1,e2,e3,e4])

q1.setDivision(nx,ny)

geoB.process()

# # Lege finite elemente aus Geometrie Quads an
fesys.getMeshCommands().getElementCommands().addFace(1,q1.getQuadList())

# # Lege Material fest
fesys.getMeshCommands().getMaterialFormulations().addMA1_2D_LinearElastic_Isotop(1,E=100,nu=0.3,thickness=1,plainstrain=0)
# #Lege FE Formulierung fest.
fesys.getMeshCommands().getElementFormulations().addEL206_Plate(num=1,mode=2,meshIdDisp=1,meshIdRot=2,dispOrder=1,rotOrder=1,EI=100,GA=50,GI=50)
# #Ersetelle Materialsatz aus FE-Formulierung und Materialformulierung1
fesys.getMeshCommands().addMaterial(1,matFormNum=1,elemFormNum=1)

fesys.getMeshCommands().setDegreesOfFreedom()

# # # Setze Randbedinungen, hier mit Hilfe von Vertices
fesys.getMeshCommands().getBoundaryConditions().singleBC(fesys.getMeshCommands().getGeometryCommands().edgeType(),number=e4.EdgeList,meshId=1,dofs=[1,1,1],shapeOrder=1)
# # #fesys.getMeshCommands().getBoundaryConditions().singleBC(fesys.getMeshCommands().getGeometryCommands().edgeType(),number=e1.EdgeList,meshId=2,dofs=[1,1,1],shapeOrder=1,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.HDiv)

# # fesys.getMeshCommands().getBoundaryConditions().singleLoad(fesys.getMeshCommands().getGeometryCommands().edgeType(),e2.EdgeList,meshId=meshid,load=[0,1,0],propnum=0,shapeorder=order)
fesys.getMeshCommands().getBoundaryConditions().singleLoad(fesys.getMeshCommands().getGeometryCommands().edgeType(),number=e2.EdgeList,meshId=1,load=[0,0,1],propnum=0,shapeorder=1,shapeType=HierAMuSPyFEM.Geometry.ShapeFunctionTypes.H1)


fesys.getMacroCommands().sparseSetUp()

# # Lastfunktion zur Skalierung der Lasten. In Abhaengigkeit der Zeit, Inkremente ueber Zeitinkrement dt. (hier "Pseudozeit")
# fesys.getMacroCommands().setPropFunction(number=0)
# fesys.getMacroCommands().setDt(1)
# fesys.getMacroCommands().timeincr()

# # # Entspricht tangent() - Befehl
# fesys.getMacroCommands().assembleSolve()

# sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(),geomNumber=v3.geoId,meshId=1)
# print(sol)
# sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(),geomNumber=v3.geoId,meshId=1)
# print(sol)

# fesys.getMacroCommands().computeEigenValues(12,24)

fesys.getMacroCommands().setLogLevel(fesys.FullLog(), fesys.FullLog())
fesys.getMacroCommands().printInfo()
# # # Entspricht  plot(paraview())
# fesys.getPlotCommands().toFile()

# # fesys.getMacroCommands().printSpMat()
