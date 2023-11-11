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
Lk.append(sympy.Float(1))
LkDiff.append(sympy.Float(0))
Lk.append(x)
LkDiff.append(sympy.Float(1))

for i in range(2,nshapes):
    faca = sympy.Rational(2*i-1,i)
    facb = sympy.Rational(i-1,i)
    temp = faca*Lk[i-1]*x - facb*Lk[i-2]
    temp = simplify(temp.evalf(dd))
    #temp = sympy.factor(temp,deep=True)
    Lk.append(temp)
    LkDiff.append(simplify(sympy.diff(temp,x)))



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

struct LegendreShapeValues{
  prec shapeValue;
  prec shapeDerivative;
};

class LegendreShapes {
public:
  LegendreShapes();
  ~LegendreShapes();
  static void getShape(prec &shape, prec &shapeDeriv, prec coor, indexType functionNumber);
  static auto getShape(prec coor, indexType functionNumber) -> LegendreShapeValues; 


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

#include <shapefunctions/LegendreShapes.h>
#include <fstream>
#include <math.h>

namespace HierAMuS {

LegendreShapes::LegendreShapes() {
  this->init = false;
  this->maxFuntions = 0;
}


LegendreShapes::~LegendreShapes() {

}


void LegendreShapes::getShape(prec &shape, prec &shapeDeriv, prec coor, indexType functionNumber) {

  switch(functionNumber) {
"""

print(type(x*2+2) == sympy.Mul)

for i in range(len(Lk)):
    a = SympyFunctionToCpp()
    a.parseFunction(Lk[i])

    temp = "  case " + str(i) + ":{\n"
    temp += "{\n" + a.getCppCode2("shape","shapeDeriv")  + "\n}\n"
    # temp += "{\n" + a.getCppCode("shape") + "\n}\n"
    # a.parseFunction(LkDiff[i])
    # temp += "{\n" + a.getCppCode("shapeDeriv") + "\n}\n"

    temp += "    }break;"

    cpptxt += temp + "\n"

cpptxt += """
default:
    throw std::runtime_error("LegendreShapes::getShape: functionNumber out of range");
}

}

auto LegendreShapes::getShape(prec coor, indexType functionNumber) -> LegendreShapeValues {

  LegendreShapeValues shape;
  switch(functionNumber) {
"""

for i in range(len(Lk)):
    a = SympyFunctionToCpp()
    a.parseFunction(Lk[i])

    temp = "  case " + str(i) + ":{\n"
    temp += "{\n" + a.getCppCode2("shape.shapeValue","shape.shapeDerivative") + "\n}\n"
    # temp += "{\n" + a.getCppCode("shape.shapeValue") + "\n}\n"
    # a.parseFunction(LkDiff[i])
    # temp += "{\n" + a.getCppCode("shape.shapeDerivative") + "\n}\n"

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

headername = pathname + "/LegendreShapes.h"
cppname = pathname + "/LegendreShapes.cpp"

hf=open(headername, 'w')
hf.write(headertxt)
hf.close()

cf=open(cppname, 'w')
cf.write(cpptxt)
cf.close()