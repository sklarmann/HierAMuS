# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import sympy


def parserSym(expr,depth=0):
    temp = ""
    depth += 1
    if len(expr.args) == 0:
        if type(expr) == sympy.Float:
            #temp += "prec(" + str(expr) + ")"
            temp += str(expr)
        else:
            temp += str(expr)
        
    else:
        if expr.func == sympy.Pow:
            t1 = parserSym(expr.args[0])
            t2 = parserSym(expr.args[1])
            temp += "pow(" + t1 + "," + t2 + ")"
        else:
            tt = []
            for i in expr.args:
                tt.append(parserSym(i,depth))

            symbol = ""
            if expr.func == sympy.Add:
                symbol = "+"
            elif expr.func == sympy.Mul:
                symbol = "*"

            if len(tt) > 1 and symbol == "+" and depth > 1:
                temp += "("
            temp += tt[0]

            

            for i in tt[1:]:
                tsym = symbol
                if tsym == "+" and i[0] == "-":
                    tsym = "-"
                temp += tsym + i.lstrip("-")
            
            if len(tt) > 1 and symbol == "+" and depth > 1:
                temp += ")"
    return temp



def valsToCppCode(vals,var):
    cpptxt = ""
    prevexp = 0
    varnames = []
    for i in range(len(vals)):
        temp = "t" + str(i)
        dval = vals[i][0]-prevexp
        prevexp = vals[i][0]
        varnames.append([temp,dval])

    temp = ""
    if varnames[0][1] != 0:
        temp = "prec " + varnames[0][0] + " = " + var + ";"
    else:
        temp = "prec " + varnames[0][0] + " = 1" + ";"
    

    for i in range(1,len(varnames)):
        temp = "prec " + varnames[i][0] + " = " + varnames[i-1][0]
        for j in range(varnames[i][1]):
            temp += "*" + var
        temp += ";"
        cpptxt += temp + "\n"


    temp = "val = "
    for i in range(len(varnames)):
        fac = str(vals[i][1]) + "*" + varnames[i][0]
        if vals[i][1] < 0:
            temp += fac
        else:
            temp += "+" + fac

    temp += ";\n\n"
    temp+= "return val;\n"
    temp+="}\n"    
    cpptxt+= temp


def symparser(texpr):
    expr = sympy.expand(texpr)

    valexp = []
    varname = ""

    if len(expr.args) == 0:
        if type(expr) != sympy.Symbol:
            valexp.append((0,expr))
        else:
            valexp.append((1,1))
    for i in expr.args:
        if len(i.args) == 0:
            valexp.append((0,i))
        else:
            if len(i.args[1].args) == 0:
                valexp.append((1,i.args[0]))
            else:
                valexp.append((i.args[1].args[1],i.args[0]))
                varname = str(i.args[1].args[0])

    
    valexp.sort(key=lambda x: x[0])

    print(valexp)
    valsToCppCode(valexp,varname)
    return valexp


class SympyFunctionToCpp:
    def __init__(self) -> None:
        self.valexp = [] # contains exponent
        self.varname = ""

        self.multvariable = [] # contains variable name and expression
        self.cppParts = [] # contains factor, multvariable
        pass

    def parseFunction(self,texpr):
        self.valexp = [] # contains exponent
        self.varname = ""

        self.multvariable = [] # contains variable name and expression
        self.cppParts = [] # contains factor, multvariable


        expr = sympy.expand(texpr)
        varname = ""

        if len(expr.args) == 0:
            if type(expr) != sympy.Symbol:
                self.valexp.append((0,expr))
            else:
                self.valexp.append((1,1))
                self.varname=str(expr)

        print(expr.func)
        if type(expr) == sympy.Mul:
            tex = expr.args[1]
            if len(tex.args) == 0:
                self.valexp.append((1,expr.args[0]))
                self.varname = str(expr.args[1])
            else:
                self.valexp.append((tex.args[1],expr.args[0]))
        else:
            for i in expr.args:
                if len(i.args) == 0:
                    self.valexp.append((0,i))
                else:
                    if len(i.args[1].args) == 0:
                        self.valexp.append((1,i.args[0]))
                        self.varname = str(i.args[1])
                    else:
                        self.valexp.append((i.args[1].args[1],i.args[0]))
                        self.varname = str(i.args[1].args[0])

        self.createVariables()

    def createVariables(self):
        self.valexp.sort(key=lambda x: x[0])
        print(self.valexp)
        pass


    def getCppCode(self,evalVarName):
        cpptxt = ""
        prevexp = 0
        varnames = []
        for i in range(len(self.valexp)):
            temp = "t" + str(i)
            dval = self.valexp[i][0]-prevexp
            prevexp = self.valexp[i][0]
            varnames.append([temp,dval])

        temp = ""
        if varnames[0][1] != 0:
            temp = "prec " + varnames[0][0] + " = " + self.varname + ";"
        else:
            temp = "prec " + varnames[0][0] + " = 1" + ";"

        cpptxt += temp + "\n"
        for i in range(1,len(varnames)):
            temp = "prec " + varnames[i][0] + " = " + varnames[i-1][0]
            for j in range(varnames[i][1]):
                temp += "*" + self.varname
            temp += ";"
            cpptxt += temp + "\n"


        temp = ""
        
        for i in range(len(varnames)):
            fac = str(self.valexp[i][1]) + "*" + varnames[i][0]
            if self.valexp[i][1] < 0:
                temp += fac
            else:
                if not temp:
                    temp = fac
                else:
                    temp += "+" + fac

         
        cpptxt+= evalVarName + " = " + temp +";"

        return cpptxt

    def getCppCode2(self,varname, diffvarname):
        cpptxt = ""
        prevexp = 0
        varnames = []
        maxExp = self.valexp[len(self.valexp)-1][0]

        expvariables = []

        for i in range(maxExp+1):
            temp = "t" + str(i)
            expvariables.append([i,temp])

        temp = ""
        temp += "prec " + expvariables[0][1] + " = 1;\n"
        for i in range(1,maxExp+1):
            temp += "prec " + expvariables[i][1] + " = " + expvariables[i-1][1] + "*" + self.varname + ";" + "\n"

        cpptxt += temp

        temp = ""
        tempdiff = ""

        for i in range(len(self.valexp)):
            exponent = self.valexp[i][0]
            if not temp:
                temp = str(self.valexp[i][1]) + "*" + expvariables[self.valexp[i][0]][1]

                
                if exponent == 0:
                    tempdiff = ""
                else:
                    temp = str(self.valexp[i][1]) + "*" + expvariables[self.valexp[i][0]][1]
                    tempdiff = str(self.valexp[i][1]*exponent) + "*" + expvariables[exponent-1][1]
            else:
                if self.valexp[i][1] < 0:
                    temp += str(self.valexp[i][1]) + "*" + expvariables[self.valexp[i][0]][1]
                    tempdiff += str(self.valexp[i][1]*exponent) + "*" + expvariables[exponent-1][1]
                else:
                    temp += "+" + str(self.valexp[i][1]) + "*" + expvariables[self.valexp[i][0]][1]
                    if not tempdiff:
                        tempdiff = str(self.valexp[i][1]*exponent) + "*" + expvariables[exponent-1][1]
                    else:
                        tempdiff += "+" + str(self.valexp[i][1]*exponent) + "*" + expvariables[exponent-1][1]
         

        cpptxt += varname + " = " + temp +";\n"
        if not tempdiff:
            tempdiff = "0"
        cpptxt += diffvarname + " = " + tempdiff +";"

        return cpptxt
