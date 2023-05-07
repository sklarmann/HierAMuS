# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from HierAMuS.HierAMuSPyWrapper import HierAMuSPyFEM

class elementFormulation:
    def __init__(self,meshcmd):
        self.mesh = meshcmd
        self.ptr = self.mesh.program.ptr
        self.elemformlist = self.ptr.getElementFormulationList()
        pass


    def addElementFormulation(self,num,element):
        self.elemformlist.addElementFormulation(num,element)

    def addEL102_Timoshenko2D(self,num,EA,EI,GA,RhoA,meshIdDisp,meshIdRot,dispOrder,rotOrder,mode,qx=0,qy=0,mz=0,propnum=0,localload=0):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("ea",EA)
        paramList.add("ei",EI)
        paramList.add("ga",GA)
        paramList.add("rhoa",RhoA)
        paramList.add("meshiddisp",meshIdDisp)
        paramList.add("meshidrot",meshIdRot)
        paramList.add("disporder",dispOrder)
        paramList.add("rotorder",rotOrder)
        paramList.add("mode",mode)
        paramList.add("qx",qx)
        paramList.add("qy",qy)
        paramList.add("mz",mz)
        paramList.add("propnum",propnum)
        paramList.add("local",localload)

        self.elemformlist.addElementFormulation(self.mesh.program.ptr,num,102,paramList)
        
    def addEL103_Timoshenko3D(self,num,E=0,G=0,A=0,Ix=0,Iy=0,Iz=0,ky=0,kz=0,geo=0,thermal=0,fe2=0,alpha=0,kappa=0,meshidDisp=1,meshidRot=2,meshidTemp=3,dispOrder=1,rotOrder=1,tempOrder=1,RVECmat=0):
        """3D Timoshenko beam as a line element. 

        Args:
            num (int): Element formulation number used to assign formulation to material set
            E (float): Young's modulus [N/m²]
            G (float): shear modulus [N/m²]
            A (float): Cross-section area [m²]
            Ix (float): St. Venant moment of area about x-axis [m^4]
            Iy (float): second moment of area about y-axis [m^4]
            Iz (float): second moment of area about z-axis [m^4]
            ky (float): shear correction factor y-direction [-]
            kz (float): shear correction factor z-direction [-]
            geo (int): 0: geometric linear, 1: finite rotations
            thermal(int): 0: basic linear elastic, 1: thermomechanically coupled (steady state)
            fe2 (int): 0: Using given parameters, 1: using FE2 scheme, specified in the material
            alpha (float): Thermal expansion coefficient [m/(m °C)]
            kappa (float): Thermal diffusity parameter [m²/s]
            meshidDisp (int): Id for displacement degrees of freedom in the mesh
            meshidRot (int): Id for rotation degrees of freedom in the mesh
            meshidTemp (int): Id for thermal degrees of freedom in the mesh
            dispOrder (int): Order of approximation of the displacement field
            rotOrder (int): Order of approximation of the rotation field
            tempOrder (int): Order of approximation of the temperature field
            RVECmat (int): If 1, uses an RVE to compute linear elastic material properties, in this case fe2 must be 0
        """
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("E",E)
        paramList.add("G",G)
        paramList.add("A",A)
        paramList.add("Ix",Ix)
        paramList.add("Iy",Iy)
        paramList.add("Iz",Iz)
        paramList.add("ky",ky)
        paramList.add("kz",kz)
        paramList.add("alpha",alpha)
        paramList.add("kappa",kappa)
        
        paramList.add("geo",geo)
        paramList.add("thermal",thermal)
        paramList.add("fe2",fe2)
        
        paramList.add("meshidDisp",meshidDisp)
        paramList.add("meshidRot",meshidRot)
        paramList.add("meshidTemp",meshidTemp)
        
        paramList.add("dispOrder",dispOrder)
        paramList.add("rotOrder",rotOrder)
        paramList.add("tempOrder",tempOrder)
        
        
        paramList.add("RVECmat",RVECmat)

        self.elemformlist.addElementFormulation(self.mesh.program.ptr,num,103,paramList)


    def addEL201_2DShell(self,num,meshiddisp,disporder,mode):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("meshiddisp",meshiddisp)
        paramList.add("disporder",disporder)
        paramList.add("mode",mode)

        self.elemformlist.addElementFormulation(self.mesh.program.ptr,num,201,paramList)

    def addEL203_BeamInterface2D(self,num,meshiddisp,meshidrot,meshidwarp,intorder,edgelist,matlist,beamVertex,mode=6):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("meshiddisp",meshiddisp)
        paramList.add("meshidrot",meshidrot)
        paramList.add("meshidwarp",meshidwarp)
        paramList.add("intorder",intorder)
        paramList.add("vertex",beamVertex)
        paramList.add("mode",mode)

        if len(matlist) != len(edgelist):
            print("Error in input of addEL203_BeamInterface2D, matlist length does not match edgelist length")
            return

        paramList.add("edgelist",edgelist)
        paramList.add("matlist",matlist)

        self.elemformlist.addElementFormulation(self.mesh.program.ptr,num,203,paramList)

    def addEL205_HDivTest(self,num,plainstrain,disporder,stressorder,mode,meshiddisp,meshidstress,E,nu):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("plainstrain",plainstrain)
        paramList.add("disporder",disporder)
        paramList.add("stressorder",stressorder)
        paramList.add("mode",mode)
        paramList.add("meshiddisp",meshiddisp)
        paramList.add("meshidstress",meshidstress)
        paramList.add("E",E)
        paramList.add("nu",nu)

        self.elemformlist.addElementFormulation(self.mesh.program.ptr,num,205,paramList)

    def addEL206_Plate(self,num,mode,meshIdDisp, meshIdRot, dispOrder, rotOrder, EI, GA, GI):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("meshiddisp",meshIdDisp)
        paramList.add("meshidrot",meshIdRot)
        paramList.add("orderdisp",dispOrder)
        paramList.add("orderrot",rotOrder)
        paramList.add("mode",mode)
        paramList.add("EI",EI)
        paramList.add("GA",GA)
        paramList.add("GI",GI)

        self.elemformlist.addElementFormulation(self.mesh.program.ptr,num,206,paramList)
    
    
    def addEL207_FaceConstraint(self,num,mode,meshIdDisp, meshIdRot, meshIdLam, meshIdMu, dispOrder, k):
        """_summary_

        Args:
            num (int): Number of the element formulation in the model
            mode (int): Currently only mode 1 is implemented (geometrically linear)
            meshIdDisp (int): Mesh Id of the displacement field
            meshIdRot (int): Mesh Id of the rotation at the point
            meshIdLam (int): Mesh Id of the translational Lagrange multiplier
            meshIdMu (int): Mesh Id of the rotational Lagrange multiplier
            dispOrder (int): Order of the displacement field
            k (float): In case of heterogeneous material, the material parameter k can be used to consider stiffness changes
        """
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("meshIdDisp",meshIdDisp)
        paramList.add("meshIdRot",meshIdRot)
        paramList.add("meshIdLam",meshIdLam)
        paramList.add("meshIdMu",meshIdMu)
        paramList.add("mode",mode)
        paramList.add("dispOrder",dispOrder)
        paramList.add("k",k)

        self.elemformlist.addElementFormulation(self.mesh.program.ptr,num,207,paramList)

    def addEL300_3DSolid(self,num,meshiddisp,disporder,mode):
        """
        Adds a 3D solid element formulation to the element formulation.

        Args:
            num: Elementformulation number used in addMaterial.
            meshiddisp: Mesh ID of the displacement field.
            disporder: Order of approximation of the displacement field.
            mode: 1: Geometrical linear, 2: Geometrical nonlinear.
        """
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("meshiddisp",meshiddisp)
        paramList.add("disporder",disporder)
        paramList.add("mode",mode)

        self.elemformlist.addElementFormulation(self.mesh.program.ptr,num,300,paramList)


    def addEL302_BeamCoupling3D(self,num,disporder,meshiddisp,meshidrot,mode=1,warptype=1,warpBounNodes=[]):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("meshiddisp",meshiddisp)
        paramList.add("meshidrot",meshidrot)
        paramList.add("disporder",disporder)
        paramList.add("mode",mode)
        paramList.add("warptype",warptype)
        
        if len(warpBounNodes) != 0:
            paramList.add("warpBounNodes",warpBounNodes)
        
        self.elemformlist.addElementFormulation(self.mesh.program.ptr,num,302,paramList)
        
    
    def addEL303_ThermoMechanikSolid3D(self,num,meshiddisp,meshidtemperatur,shapeorder,mu,lamb,alpha,c,rho0,T0,kappa,mode=1):
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("meshiddisp",meshiddisp)
        paramList.add("meshidtemperatur",meshidtemperatur)
        paramList.add("disporder",shapeorder)
        paramList.add("mu",mu)
        paramList.add("lambda",lamb)
        paramList.add("alpha",alpha)
        paramList.add("c",c)
        paramList.add("rho0",rho0)
        paramList.add("T0",T0)
        paramList.add("kappa",kappa)
        
        
        paramList.add("mode",mode)

        self.elemformlist.addElementFormulation(self.mesh.program.ptr,num,303,paramList)
        
    def addEL307_VolumeConstraint(self,num,meshiddisp,meshIdLam,shapeorder,stiffness,center,mode=1):
        """_summary_

        Args:
            num (int): Number of the element formulation in the model
            meshiddisp (int): Mesh Id of the displacement field
            meshIdLam (int): Mesh Id of the Lagrange multiplier
            shapeorder (int): Shape order of the displacement field.
            stiffness (float): Relative stiffness of the constraint, required if different materials are used.
            center (int): 0: Moment has max value at the boundary of the element; 1: Moment has max value at the center of the element.
            mode (int, optional): Mode 1: Constrain linear moment in x-Direction (Beam RVE). Defaults to 1.
        """
        paramList = HierAMuSPyFEM.ParameterList()
        paramList.add("meshiddisp",meshiddisp)
        paramList.add("disporder",shapeorder)
        paramList.add("meshIdLam",meshIdLam)
        paramList.add("stiffness",stiffness)
        paramList.add("center",center)
        
        paramList.add("mode",mode)

        self.elemformlist.addElementFormulation(self.mesh.program.ptr,num,307,paramList)
