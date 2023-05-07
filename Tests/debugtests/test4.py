# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import sympy
from sympy.simplify.simplify import simplify
from sympy.utilities.codegen import codegen

dd = 3
nshapes = 2

x=sympy.Symbol('x')


one = sympy.N(1,dd)
two = sympy.N(2,dd)
three = sympy.N(3,dd)


Lk= []
Lk.append(sympy.Rational(1,3))
Lk.append(x)

for i in range(2,nshapes):
    faca = sympy.Rational(2*i-1,i)
    facb = sympy.Rational(i-1,i)
    temp = faca*Lk[i-1]*x - facb*Lk[i-2]
    temp = simplify(temp)
    #temp = sympy.separatevars(temp,force=True)
    Lk.append(temp)


# for i in Lk:
#     temp = sympy.solve(sympy.N(i,dd),x)
#     print(temp)
#     print(sympy.N(temp,dd))
#     print(sympy.solve(sympy.N(i,dd),x,numerical=True,rational=False))



from sympy.utilities.codegen import codegen


# for i in Lk:
#     ##print(sympy.printing.ccode(i.evalf(dd),precision=5))
#     for j in i.as_coeff_Add():
#         print(j)
#     # [(c_name, c_code), (h_name, c_header)] = codegen(

#     #     ("f", (i.evalf(dd))), "C89", "test", header=False, empty=True)
#     # #print(c_name)
#     # print(i)
#     # print(i.evalf(dd))
#     # print(c_code)


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
                temp += symbol + i
            
            if len(tt) > 1 and symbol == "+" and depth > 1:
                temp += ")"
    return temp

f = (x+2)**(2)+x+sympy.Rational(1,3)
#f = sympy.simplify(f)
#print(Lk[4])
print(parserSym(f.evalf(3)))
print(f.evalf(3))

# LkDiff = []
# for i in Lk:
#     temp = i.diff(x)
#     LkDiff.append(temp)

# for i in range(len(Lk)):
#     print("Legendre"+str(i)+"(prec x){")
#     print("\treturn " + parserSym((Lk[i]).evalf(dd)) + ";")
#     print("}\n")

# for i in range(len(Lk)):
#     print("LegendreDiff"+str(i)+"(prec x){")
#     print("\treturn " + parserSym((LkDiff[i]).evalf(dd)) + ";")
#     print("}\n")