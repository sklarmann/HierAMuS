# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import cppyy

class OutputHandler:
    def __init__(self):
        self.LogFile = ''
        self.currentPrintLevel = 0
        self.currentWriteLevel = 0
        pass

    def setLogFile(self,fileName):
        self.LogFile = fileName
        self.File = open(self.LogFile,'w')

    def setLogLevel(self,printIn,writeIn):
        self.currentPrintLevel=printIn
        self.currentWriteLevel=writeIn

    def output(self,printLevel,writeLevel,*args):
        if printLevel <= self.currentPrintLevel or writeLevel <= self.currentWriteLevel:
            ostr = ''
            for i in args:
                ostr += str(i)
            
            if printLevel <= self.currentPrintLevel:
                print(ostr)

            if writeLevel <= self.currentWriteLevel:
                pass




a=OutputHandler()
a.setLogLevel(1,1)
b=cppyy.gbl.std.vector['int']()
b.push_back(1)
a.output(0,0,"blabla took:  ",2, " seconds", b)