#!/usr/bin/env python

# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/klarmann/git/cppfemtemp')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Box_1 = geompy.MakeBoxDXDYDZ(200, 200, 200)
[Face_1,Face_2,Face_3,Face_4,Face_5,Face_6] = geompy.ExtractShapes(Box_1, geompy.ShapeType["FACE"], True)
[Face_1, Face_2, Face_3, Face_4, Face_5, Face_6] = geompy.GetExistingSubObjects(Box_1, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Box_1, 'Box_1' )
geompy.addToStudyInFather( Box_1, Face_1, 'Face_1' )
geompy.addToStudyInFather( Box_1, Face_2, 'Face_2' )
geompy.addToStudyInFather( Box_1, Face_3, 'Face_3' )
geompy.addToStudyInFather( Box_1, Face_4, 'Face_4' )
geompy.addToStudyInFather( Box_1, Face_5, 'Face_5' )
geompy.addToStudyInFather( Box_1, Face_6, 'Face_6' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Box_1)
Regular_1D = Mesh_1.Segment()
Number_of_Segments_1 = Regular_1D.NumberOfSegments(6)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Hexa_3D = Mesh_1.Hexahedron(algo=smeshBuilder.Hexa)
Face_1_1 = Mesh_1.GroupOnGeom(Face_1,'Face_1',SMESH.FACE)
Face_2_1 = Mesh_1.GroupOnGeom(Face_2,'Face_2',SMESH.FACE)
Face_3_1 = Mesh_1.GroupOnGeom(Face_3,'Face_3',SMESH.FACE)
Face_4_1 = Mesh_1.GroupOnGeom(Face_4,'Face_4',SMESH.FACE)
Face_5_1 = Mesh_1.GroupOnGeom(Face_5,'Face_5',SMESH.FACE)
Face_6_1 = Mesh_1.GroupOnGeom(Face_6,'Face_6',SMESH.FACE)
isDone = Mesh_1.Compute()
[ Face_1_1, Face_2_1, Face_3_1, Face_4_1, Face_5_1, Face_6_1 ] = Mesh_1.GetGroups()
Volumes = Mesh_1.CreateEmptyGroup( SMESH.VOLUME, 'Volumes' )
Volumes.AddFrom( Mesh_1.GetMesh() )
loadedVerts = Mesh_1.GroupOnGeom(Face_6,'LoadVerts',SMESH.NODE)

AllVerts = Mesh_1.CreateEmptyGroup( SMESH.NODE, 'AllVerts' )
AllVerts.AddFrom( Mesh_1.GetMesh() )

## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Hexa_3D.GetAlgorithm(), 'Hexa_3D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Face_1_1, 'Face_1')
smesh.SetName(Face_2_1, 'Face_2')
smesh.SetName(Face_3_1, 'Face_3')
smesh.SetName(Face_4_1, 'Face_4')
smesh.SetName(Face_5_1, 'Face_5')
smesh.SetName(Face_6_1, 'Face_6')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()


import HierAMuS

fesys = HierAMuS.FEMPy("/home/klarmann/git/cppfemtemp", "test.log")
fesys.setStaticSolutionState()
fesys.setSolver(4)
fesys.initSalome(SMESH, smesh)
fesys.getMacroCommands().setLogLevel(fesys.FullLog(), fesys.FullLog())

fesys.getMeshCommands().getFromSalome().createInternalEdges(Mesh_1)
fesys.getMeshCommands().getFromSalome().geomFromSalome(Mesh_1)

fesys.getMeshCommands().getFromSalome().addElementFromSalomeGroup(Volumes, 1)

fesys.getMeshCommands().getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1, 100, 0.3)
fesys.getMeshCommands().getElementFormulations().addEL300_3DSolid(1, 1, 1, 2)
fesys.getMeshCommands().addMaterial(1, 1, 1)

fesys.getMeshCommands().setDegreesOfFreedom()

fesys.getMeshCommands().getBoundaryConditions().singleBC(fesys.getMeshCommands().getGeometryCommands().faceType(), number=Face_1_1.GetIDs(), meshId=1, dofs=[1,1,1], shapeOrder=1)
fesys.getMeshCommands().getBoundaryConditions().singleLoad(fesys.getMeshCommands().getGeometryCommands().faceType(), number=Face_6_1.GetIDs(), meshId=1, load=[1,0,0], propnum=0)

fesys.getMacroCommands().sparseSetUp()

fesys.getMacroCommands().setPropFunction(0)
fesys.getMacroCommands().setDt(50)
fesys.getMacroCommands().timeincr()

fesys.getMacroCommands().newton()

# for i in AllVerts.GetIDs():
#   sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(), i, 1)
#   print([i," ",sol])

a=loadedVerts.GetIDs()[0]
sol = fesys.getMacroCommands().getSolution(fesys.getMeshCommands().getGeometryCommands().vertexType(), a, 1)
print([a," ",sol])

fesys.getPlotCommands().toFile()
#fesys.getMacroCommands().computeEigenValues(4,8)
#fesys.getMacroCommands().printInfo()