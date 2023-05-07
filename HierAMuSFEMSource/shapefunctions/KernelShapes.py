# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import sympy
from sympy.simplify.simplify import simplify
from sympy.utilities.codegen import codegen


# Größtenteils identisch mit Lobattoshapes
# Kernelshapes Ker = Lob / (Lob[0]*Lob[1])


dd = 30   #Rechengenauigkeit des Skripts, ist höher da Polynomdivision genauigkeit verschlechtert
ddd = 20  #Ausgegebene Nachkommastellen
nshapes = 20

x=sympy.Symbol('coor')


one = sympy.N(1,dd)
two = sympy.N(2,dd)
three = sympy.N(3,dd)


Lk= []
LkDiff = []
Lob = []
Lk.append(sympy.Float(1))
Lk.append(x)

for i in range(2,nshapes):
    faca = sympy.Rational(2*i-1,i)
    facb = sympy.Rational(i-1,i)
    temp = faca*Lk[i-1]*x - facb*Lk[i-2]
    temp = simplify(temp.evalf(dd))
    #temp = sympy.factor(temp,deep=True)
    Lk.append(temp)

Lob.append(((1-x)/2).evalf(dd))
Lob.append(((1+x)/2).evalf(dd))
Ker = []

for i in range(2,nshapes):
  faca = sympy.Rational(2,2*i-1)
  faca = sympy.sqrt(faca).evalf(dd)
  Ker.append(sympy.simplify(sympy.integrate(Lk[i-1],(x,-1,x))/faca))
  q, r = sympy.div(Ker[i-2], Lob[0]*Lob[1])
  print(r)
  Ker[i-2] = q.evalf(ddd)

print(Ker)

Lk=[]
Lk = Ker

# for i in range(nshapes):
#   LkDiff.append(sympy.simplify(sympy.diff(Lk[i],x)))

from parser1 import SympyFunctionToCpp

import os, sys

print('sys.argv[0] =', sys.argv[0])             
pathname = os.path.dirname(sys.argv[0])        
print('path =', pathname)
print('full path =', os.path.abspath(pathname)) 





headertxt = """
//
// Created by simon on 7/4/21.
//

#pragma once
#include <vector>
#include <string>
#include <datatypes.h>



namespace HierAMuS {

struct KernelShapesValues{
  prec shapeValue;
  prec shapeDerivative;
};

class KernelShapes {
public:
  KernelShapes();
  ~KernelShapes();
  void readData(std::string &filename);
  static void getShape(prec &shape, prec &shapeDeriv, prec coor, indexType functionNumber);
  static auto getShape( prec coor, indexType functionNumber) -> KernelShapesValues;


private:
  bool init;
  indexType maxFuntions;
  std::vector<std::vector<prec>> coefficients;
  std::vector<std::vector<indexType>> exponents;

};


}
"""

cpptxt = """
//
// Created by simon on 7/4/21.
//
#include <preprocessDefine.h>
#include <shapefunctions/KernelShapes.h>
#include <fstream>
#include <math.h>

namespace HierAMuS {

KernelShapes::KernelShapes() {
  this->init = false;
  this->maxFuntions = 0;
}

KernelShapes::~KernelShapes() {

}

void KernelShapes::readData(std::string &filename) {

}

void KernelShapes::getShape(prec &shape, prec &shapeDeriv, prec coor, indexType functionNumber) {

  switch(functionNumber) {
"""


for i in range(len(Lk)):
    a = SympyFunctionToCpp()
    a.parseFunction(Lk[i])

    temp = "  case " + str(i) + ":{\n"
    temp += "{\n" + a.getCppCode2("shape","shapeDeriv")  + "\n}\n"

    temp += "    }break;"

    cpptxt += temp + "\n"

cpptxt += """
default:
    throw std::runtime_error("LegendreShapes::getShape: functionNumber out of range");
}

}

auto KernelShapes::getShape(prec coor, indexType functionNumber) -> KernelShapesValues {

  KernelShapesValues shape;
  switch(functionNumber) {
"""

for i in range(len(Lk)):
    a = SympyFunctionToCpp()
    a.parseFunction(Lk[i])

    temp = "  case " + str(i) + ":{\n"
    temp += "{\n" + a.getCppCode2("shape.shapeValue","shape.shapeDerivative") + "\n}\n"


    temp += "    }break;"

    cpptxt += temp + "\n"


cpptxt += """
default:
    throw std::runtime_error("LegendreShapes::getShape: functionNumber out of range");
}

  return shape;
}

}
"""

# print(cpptxt)

headername = pathname + "/KernelShapes.h"
cppname = pathname + "/KernelShapes.cpp"

hf=open(headername, 'w')
hf.write(headertxt)
hf.close()

cf=open(cppname, 'w')
cf.write(cpptxt)
cf.close()

# print(Lob[1])
# print(Lk[1])
# print(Lk[2])
# print(Ker[1])
# print(Ker[2])

b = Lk[0]

a = SympyFunctionToCpp()
a.parseFunction(b)

# print("break")
# print(b)
# print(sympy.diff(b,x))
# print(a.getCppCode2("shape","shapederiv"))