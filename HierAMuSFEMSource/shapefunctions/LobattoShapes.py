# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import sympy
from sympy.simplify.simplify import simplify
from sympy.utilities.codegen import codegen

dd = 20
nshapes = 30

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

for i in range(2,nshapes):
  faca = sympy.Rational(2,2*i-1)
  faca = sympy.sqrt(faca).evalf(dd)
  Lob.append(sympy.simplify(sympy.integrate(Lk[i-1],(x,-1,x))/faca))


Lk=[]
Lk = Lob

for i in range(nshapes):
  LkDiff.append(sympy.simplify(sympy.diff(Lk[i],x)))

from parser import parserSym, symparser, SympyFunctionToCpp

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

struct LobattoShapesValues{
  prec shapeValue;
  prec shapeDerivative;
};

class LobattoShapes {
public:
  LobattoShapes();
  ~LobattoShapes();
  static void getShape(prec &shape, prec &shapeDeriv, prec coor, indexType functionNumber);
  static auto getShape( prec coor, indexType functionNumber) -> LobattoShapesValues;


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

#include <shapefunctions/LobattoShapes.h>
#include <fstream>
#include <math.h>

namespace HierAMuS {

LobattoShapes::LobattoShapes() {
  this->init = false;
  this->maxFuntions = 0;
}

LobattoShapes::~LobattoShapes() {

}


void LobattoShapes::getShape(prec &shape, prec &shapeDeriv, prec coor, indexType functionNumber) {

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

auto LobattoShapes::getShape(prec coor, indexType functionNumber) -> LobattoShapesValues {

  LobattoShapesValues shape;
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

print(cpptxt)

headername = pathname + "/LobattoShapes.h"
cppname = pathname + "/LobattoShapes.cpp"

hf=open(headername, 'w')
hf.write(headertxt)
hf.close()

cf=open(cppname, 'w')
cf.write(cpptxt)
cf.close()

print(Lob[1])
print(Lob[2])
print(Lk[1])
print(Lk[2])

b = Lk[0]

a = SympyFunctionToCpp()
a.parseFunction(b)

print("break")
print(b)
print(sympy.diff(b,x))
print(a.getCppCode2("shape","shapederiv"))