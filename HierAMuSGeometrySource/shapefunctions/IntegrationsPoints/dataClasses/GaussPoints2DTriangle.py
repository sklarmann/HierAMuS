# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import sympy
from sympy.simplify.simplify import simplify
from sympy.utilities.codegen import codegen
import math


dd = 5
nshapes = 3

x=sympy.Symbol('x')
y=sympy.Symbol('y')


one = sympy.N(1,dd)
two = sympy.N(2,dd)
three = sympy.N(3,dd)

Fak = []

for i in range(0,nshapes+5):
    Fak.append(math.factorial(i))

Int = []
for n in range(1,nshapes+1):
    I = []
    for i in range(0,n):
        for j in range(0,n-i):
            I.append(Fak[i]*Fak[j]/(Fak[i+j+2]))
    Int.append(I)
print(Int)

