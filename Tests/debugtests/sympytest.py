# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import sympy
from sympy.simplify.simplify import simplify
from sympy.utilities.codegen import codegen

N1x = sympy.Symbol('N1x')
N1y = sympy.Symbol('N1y')
N1z = sympy.Symbol('N1z')

N2x = sympy.Symbol('N2x')
N2y = sympy.Symbol('N2y')
N2z = sympy.Symbol('N2z')

B1 = sympy.Matrix([[N1x,0,0],[0,N1y,0],[0,0,N1z],[N1y,N1x,0],[N1z,0,N1y],[0,N1z,N1x]])
B2 = sympy.Matrix([[N2x,0,0],[0,N2y,0],[0,0,N2z],[N2y,N2x,0],[N2z,0,N2y],[0,N2z,N2x]])

C=sympy.Matrix(sympy.MatrixSymbol("C",6,6))

for i in range(6):
    for j in range(6):
        C[i,j] = sympy.Symbol('C'+str(i)+str(j))

for i in range(6):
    for j in range(i+1):
        C[i,j] = C[j,i]

s = B1.transpose()*C*B2
s = sympy.simplify(s)

for i in range(3):
    for j in range(3):
        print((s[i,j]))