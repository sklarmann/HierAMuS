# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from HierAMuS.HierAMuSPyWrapper import HierAMuSPyFEM

class constraintCommands:
    def __init__(self,meshcmd):
        self.mesh = meshcmd
        self.ptr = self.mesh.program.ptr
        
        pass

    def generalLink(self,geoType,masterNumber,slaveNumber,meshId,order,masterDof,slaveDof,factor,difference,reorient=True):
        solState = self.ptr.getSolutionState()
        constHandler = solState.getConstraintHandler()
        
        if type(masterNumber) != list:
            masterNumber = [masterNumber]
            
        if type(slaveNumber) != list:
            slaveNumber = [slaveNumber]
        
        constHandler.GeneralLinkGeo(self.ptr,geoType,masterNumber,slaveNumber,meshId,order,masterDof,slaveDof,factor,difference,reorient)