# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

from CPPFEMPython import HierAMuSPyFEM


class pyPointerCollection(HierAMuSPyFEM.PointerCollection):
    def __init__(self):
        self.solution = 0
        super().__init__()

    def setSolutionState(self,solstate):
        self.solution = solstate
        self.solution.setProps()

    def getSolutionState(self):
        return self.solution


    