# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

class plotCommands:
    def __init__(self,program) -> None:
        self.prog = program
        pass

    def toFile(self):
        self.prog.ptr.getPlotControlInterface().toFile(self.prog.ptr)
        #self.prog.pyPlot.writeFile()
        #plotInterface = self.prog.ptr.getVtkPlotInterface()
        #/plotInterface = plotInterface
        #plotInterface.writeFile(pathName,fileName)