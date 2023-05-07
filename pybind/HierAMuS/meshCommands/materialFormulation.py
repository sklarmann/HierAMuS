# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from HierAMuS.HierAMuSPyWrapper import HierAMuSPyFEM

class materialFormulation:
    def __init__(self,meshcmd):
        self.mesh = meshcmd
        self.ptr = self.mesh.program.ptr
        self.matlist = self.ptr.getMaterialFormulationList()
        pass

    def addMaterialFormulation(self,number,materialForm):
        self.matlist.addMaterial(number,materialForm)

    def addMA1_2D_PlainStrain_3D(self,number, number3D):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("number",number3D)

        self.matlist.addMaterial(self.ptr,number,201,paramList)


    def addMA3_2D_LinearElastic_Isotrop(self,number,E,nu,thickness,plainstrain):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("emodul",E)
        paramList.add("nu",nu)
        paramList.add("thickness",thickness)
        paramList.add("plainstrain",plainstrain)

        self.matlist.addMaterial(self.ptr,number,203,paramList)



    def addMA1_3D_LinearElastic_Isotrop(self,number,E,nu):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("emodul",E)
        paramList.add("nu",nu)

        self.matlist.addMaterial(self.ptr,number,301,paramList)

    def addMA2_3D_NeoHook(self,number,Lambda,G):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("Lambda",Lambda)
        paramList.add("G",G)

        self.matlist.addMaterial(self.ptr,number,302,paramList)
        
        
    def addMA3_SmallStrainPlasticity(self,number,E,nu,y0,yinf,xh,xd,eta,maxiterations=20):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("emodul",E)
        paramList.add("nu",nu)
        paramList.add("y0",y0)
        paramList.add("yinf",yinf)
        paramList.add("xh",xh)
        paramList.add("xd",xd)
        paramList.add("eta",eta)
        paramList.add("maxiterations",maxiterations)

        self.matlist.addMaterial(self.ptr,number,303,paramList)

    def addMAS1_Homogenization(self,number,RVE,maxIterations=10):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("maxIterations",maxIterations)
        
        self.matlist.addMaterial(self.ptr,number,400,paramList)
        mat = self.matlist.getMaterial(number)
        mat.setRVE(self.ptr,RVE.ptr)
        
        