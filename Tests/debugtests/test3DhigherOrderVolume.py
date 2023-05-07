# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import HierAMuS

import os
import sys
import time
import gmsh
import math
import matplotlib.pyplot as plt

# time.sleep(20)



def calc(refinements,order):
    x=[]
    y=[]

    for ll in range(refinements):
        gmsh.initialize()

        gmsh.model.add("t11")

        L=10
        b=1
        h=1
        angle = 0

        transmesh = True
        fac = 2**ll
        nx = 10*fac
        ny = 1*fac

        box = gmsh.model.occ.addBox(0, 0, 0, L, b, h)
        p1 = gmsh.model.occ.addPoint(10.01, 0.5, 0.5)



        #gmsh.model.occ.rotate([(3,box),(0,p1)], 0, 0.5, 0.5, 1, 0, 0, angle/180.0*math.pi)
        print(p1)
        surfLoopTag, surfTags = gmsh.model.occ.getSurfaceLoops(box)
        back = surfTags[0][0]
        front = surfTags[0][1]


        gmsh.model.occ.synchronize()


        if transmesh:
            # Getting curves from front and back surfaces
            # Marking them as transfinite curves with ny nodes
            curveLoopTag, curveTag = gmsh.model.occ.getCurveLoops(back)
            for i in curveTag[0]:
                gmsh.model.mesh.setTransfiniteCurve(i, ny)

            curveLoopTag, curveTag = gmsh.model.occ.getCurveLoops(front)
            for i in curveTag[0]:
                gmsh.model.mesh.setTransfiniteCurve(i, ny)

            # Getting curves in length direction
            edgelists = []
            for face in surfTags[0][2:]:
                cltag, edges = gmsh.model.occ.getCurveLoops(face)
                edgelists.append(edges[0].tolist())

            a = []
            for i in range(4):
                for j in range(i+1,4):
                    a.append(set(edgelists[i]).intersection(edgelists[j]))

            for i in a:
                for j in i:
                    gmsh.model.mesh.setTransfiniteCurve(j, nx)

            for i in surfTags[0]:
                gmsh.model.mesh.setTransfiniteSurface(i)

            gmsh.model.mesh.setTransfiniteVolume(box)

        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
        gmsh.option.setNumber("Mesh.RecombineAll", 2)
        #gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2)
        gmsh.model.mesh.generate(0)
        gmsh.model.mesh.generate(1)
        gmsh.model.mesh.generate(2)
        gmsh.model.mesh.generate(3)
        #gmsh.model.mesh.refine()


        #gmsh.fltk().run()

        # Ermittle Pfad des Skripts
        print('sys.argv[0] =', sys.argv[0])
        pathname = os.path.dirname(sys.argv[0])
        print('path =', pathname)
        print('full path =', os.path.abspath(pathname))

        currPath = os.path.abspath(pathname)

        fesys = HierAMuS.FEMPy(currPath, "test3DBeamInterface")
        fesys.setStaticSolutionState()
        fesys.setSolver(4)
        #fesys.setMaxThreads(2)
        Macro = fesys.getMacroCommands()
        Macro.setLogLevel(fesys.FullLog(), fesys.NoLog())


        Mesh = fesys.getMeshCommands()
        Gm = Mesh.getFromGMESH()

        Gm.addVertices(gmsh)
        Gm.addBrickFiniteElements(gmsh, order=1, tag=box, material=1)

        Geo = Mesh.getGeometryCommands()

        # Adding the Volume Elements and assign the material
        Mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(1, E=1000, nu=0 )
        Mesh.getElementFormulations().addEL300_3DSolid(1, meshiddisp=1, disporder=order, mode=1)
        Mesh.addMaterial(1, 1, 1)

        # Adding the interface element on the load side and assign the material
        faces = Gm.getFaceNumbers(gmsh, tag=front, ftype=4, order=1)
        vnum = Geo.requestVertex(11, 0.0, 0.0)

        nt,cn,t =gmsh.model.mesh.getNodes(dim=0,tag=p1)
        #vnum = Geo.requestVertex(cn[0], cn[1], cn[2])
        vnum = nt[0]

        print(cn[0])
        #Mesh.getElementCommands().add3DInterfaceElement(materialNumber=2, specialNumber=1, vertNum=vnum, faceList=faces)
        #Mesh.getElementFormulations().addEL302_BeamCoupling3D(num=2, disporder=1, meshiddisp=1, meshidrot=2)
        #Mesh.addMaterial(2, 1, 2)


        Macro.setLogLevel(fesys.FullLog(),fesys.BasicLog())

        Mesh.setDegreesOfFreedom()

        fnums = Gm.getFaceNumbers(gmsh, tag=back, ftype=4, order=1)
        Mesh.getBoundaryConditions().singleBC(Geo.faceType(), number=fnums, meshId=1, dofs=[1,1,1], shapeOrder=order)
        fnums = Gm.getFaceNumbers(gmsh, tag=front, ftype=4, order=1)
        Mesh.getBoundaryConditions().singleLoad(Geo.faceType(), number=fnums, meshId=1, load=[0,1,0], propnum=1,shapeorder=order)
        #Mesh.getBoundaryConditions().singleLoad(Geo.vertexType(), number=vnum, meshId=2, load=[1,0,0], propnum=1)

        Macro.sparseSetUp()

        Macro.setPropFunction(1)
        Macro.setPropFunction(0)
        Macro.setDt(1)
        Macro.timeincr()


        Macro.newton(maxIteration=2)
        Macro.timeincr()
        # Macro.printInfo()

        # fesys.getPlotCommands().toFile()



        gmsh.finalize()

        sol=Macro.getSolution(Geo.vertexType(),7,1)
        x.append(ll+1)
        y.append(sol[1])
        #print(sol)
        #sol=Macro.getSolution(Geo.vertexType(),vnum,2)
        #print(sol)
    return x,y


def calcErr(ref,sol):
    err=[]
    for i in sol:
        err.append(abs((ref-i)/ref))
    ferr = err[0]
    #for i in range(len(err)):
    #    err[i]/=ferr
        
    return err

fig, ax = plt.subplots()
ref = 4.02436100177911

numref=5

x,y = calc(numref+1,1)
x2,y2 = calc(numref,2)
x3,y3 = calc(numref-1,3)
x4,y4 = calc(numref-2,4)


err1 = calcErr(ref,y)
err2 = calcErr(ref,y2)
err3 = calcErr(ref,y3)
err4 = calcErr(ref,y4)

    

ax.loglog(x, err1, linewidth=2.0)
ax.loglog(x2, err2, linewidth=2.0)
ax.loglog(x3, err3, linewidth=2.0)
ax.loglog(x4, err4, linewidth=2.0)
plt.show()
