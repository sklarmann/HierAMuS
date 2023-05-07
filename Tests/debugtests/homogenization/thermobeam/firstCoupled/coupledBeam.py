# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import os,sys
import HierAMuS
from matplotlib import pyplot as plt
import gmsh
from RVE import rectAngularRVE

RVE = rectAngularRVE(L=1,b=1,h=1,nx=1,ny=2,nz=2,order=2)
#RVE.getMacroCommands().setLogLevel(RVE.BasicLog(),RVE.BasicLog())
RVE.getMacroCommands().setLogLevel(RVE.NoLog(),RVE.NoLog())
#RVE.getMacroCommands().getHomogenizationCommands().setStrains([0.1,0.0,0.0,0.0,0.0,0.0])
#RVE.getMacroCommands().newton()
#RVE.getMacroCommands().assembleSolve()
#RVE.getMacroCommands().getHomogenizationCommands().homogenize()
#
#exit()

# current path of file
currPath = os.path.dirname(os.path.realpath(__file__))

# geometry parameters
L=10
E=100
nu=0.0
G=E/(2*(1+nu))
A=1
Ix=0.3141
Iy=1.0/12.0
Iz=Iy
ky=5.0/6.0
kz=ky
alpha=0.1
kappa=1

# mesh
nx=20
dispOrder = 3
rotOrder = 2
tempOrder = 1
meshIdDisp = 1
meshIdRot = 2
meshIdTemp = 3



fesys = HierAMuS.FEMPy(currPath,"coupledBeam")
fesys.getMacroCommands().setLogLevel(fesys.BasicLog(),fesys.BasicLog())
fesys.setStaticSolutionState()
fesys.setSolver(2)

geo = fesys.getMeshCommands().getGeometryCommands()
elem = fesys.getMeshCommands().getElementCommands()

dx = L/nx
geo.addVertex(1,0,0,0)
for i in range(nx):
    nn = i+2
    geo.addVertex(nn,(i+1)*dx,0,0)
    geo.addLinearEdge(i+1,[i+1,nn])
    elem.addEdge(1,i+1)
    
    
    
mesh = fesys.getMeshCommands()

mesh.getElementFormulations().addEL103_Timoshenko3D(num=1,E=E,G=G,A=A,Ix=Ix,Iy=Iy,Iz=Iz,ky=ky,kz=kz,geo=0,thermal=0,fe2=1,alpha=alpha,kappa=kappa,meshidDisp=meshIdDisp,meshidRot=meshIdRot,meshidTemp=meshIdTemp,dispOrder=dispOrder,rotOrder=rotOrder,tempOrder=tempOrder)
mesh.getMaterialFormulations().addMA1_3D_LinearElastic_Isotrop(number=1,E=E,nu=nu)
mesh.getMaterialFormulations().addMAS1_Homogenization(number=2,RVE=RVE)
mesh.addMaterial(matNum=1,matFormNum=2,elemFormNum=1)



mesh.setDegreesOfFreedom()

boun = mesh.getBoundaryConditions()
boun.BCVertex(number=1,meshId=meshIdDisp,dofs=[1,1,1])
boun.BCVertex(number=1,meshId=meshIdRot,dofs=[1,1,1])

boun.LoadVertex(number=nx+1,meshId=meshIdDisp,load=[0,0.05,0],propnum=1)

fesys.getMacroCommands().sparseSetUp()

def timeFunc(t):
    if t<5:
        return t
    else:
        return 10-t 

fesys.getMacroCommands().setPropFunction(1,function=timeFunc)
fesys.getMacroCommands().setDt(1)


x=[0]
y=[0]
for i in range(10):
    fesys.getMacroCommands().timeincr()

    fesys.getMacroCommands().newton(refResidual=1e-10)
    sol = fesys.getMacroCommands().getSolution(geo.vertexType(),nx+1,1)
    x.append(timeFunc(i+1))
    y.append(sol[1])
    


fig, ax = plt.subplots()
ax.plot(y, x, linewidth=2.0)
plt.show()
#fesys.getMacroCommands().assembleSolve()

print(fesys.getMacroCommands().getSolution(geo.vertexType(),nx+1,1))
    
#fesys.getMacroCommands().printInfo()