# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import sympy
from sympy.printing import ccode

# create J-Matrix

def getMatrixRow(exp,ematVoigt):
    expanded = sympy.expand(exp)

    coeffs = []
    for i in range(6):
        val = ematVoigt[i]
        coeffs.append(expanded.coeff(val,1))

    return coeffs



JField = []
for i in range(3):
    row = []
    for j in range(3):
        name = "J_"+str(i+1) + str(j+1)
        row.append(sympy.symbols(name,real=True))
    JField.append(row)

J = sympy.Matrix(JField)

nn = ["x","y","z"]
EField=[]
for i in range(3):
    row=[]
    for j in range(3):
        name = "E_" + nn[i] + nn[j]
        row.append(sympy.symbols(name,real=True))
    EField.append(row)

E = sympy.Matrix(EField)

EFieldVoigt =[]
for i in range(3):
    EFieldVoigt.append(E[i,i])
EFieldVoigt.append(E[0,1])
EFieldVoigt.append(E[0,2])
EFieldVoigt.append(E[1,2])

EFieldVoigt = sympy.Matrix(EFieldVoigt)

Emapped = J*E*sympy.transpose(J)
Emapped = sympy.simplify(Emapped)

EVoigt = []
for i in range(3):
    EVoigt.append(Emapped[i,i])

EVoigt.append(Emapped[0,1])
EVoigt.append(Emapped[0,2])
EVoigt.append(Emapped[1,2])

EVoigt = sympy.Matrix(EVoigt)

transformMat = []
for i in range(6):
    #print(getMatrixRow(EVoigt[i],EFieldVoigt))
    transformMat.append(getMatrixRow(EVoigt[i],EFieldVoigt))

J = sympy.MatrixSymbol('J0',3,3)
Res = sympy.transpose(J)*E*J
EVoigt = []
for i in range(3):
    EVoigt.append(Res[i,i])

EVoigt.append((Res[0,1]+Res[1,0]))
EVoigt.append((Res[0,2]+Res[2,0]))
EVoigt.append((Res[1,2]+Res[2,1]))

transmat = []
for i in range(6):
    transmat.append(getMatrixRow(EVoigt[i],EFieldVoigt))
    #print(getMatrixRow(EVoigt[i],EFieldVoigt))
    #print(ccode(getMatrixRow(EVoigt[i],EFieldVoigt)))

for i in range(len(transmat)):
    h = transmat[i]
    for j in range(len(h)):
        exp = h[j]
        exp = sympy.expand(exp)
        expstr = str(exp)
        expstr = expstr.replace('[','(')
        expstr = expstr.replace(']',')')
        print('trans(' + str(i) + ',' + str(j) + ') = ' + expstr + ";")