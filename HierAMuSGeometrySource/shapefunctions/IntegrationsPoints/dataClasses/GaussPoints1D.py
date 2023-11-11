# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import sympy
from sympy.simplify.simplify import simplify
from sympy.utilities.codegen import codegen
from mpmath import *


dd = 20
mp.dps = dd; mp.pretty = False
nshapes = 40

x=sympy.Symbol('coor',real=True)


one = sympy.N(1,dd)
two = sympy.N(2,dd)
three = sympy.N(3,dd)


Lk= []
LkDiff = []
Lk.append(sympy.Rational(1,1))
LkDiff.append(sympy.Rational(0,1))
Lk.append(x)
LkDiff.append(sympy.Rational(1,1))

for i in range(2,nshapes):
    faca = sympy.Rational(2*i-1,i)
    facb = sympy.Rational(i-1,i)
    temp = faca*Lk[i-1]*x - facb*Lk[i-2]
    temp = simplify(temp)
    #temp = sympy.factor(temp,deep=True)
    Lk.append(temp)
    LkDiff.append(sympy.simplify(sympy.diff(Lk[i],x)))


GP = []
Weights = []
for i in range(1,nshapes):
    t = sympy.Poly(Lk[i], x)
    coe = t.all_coeffs()
    coe2 = []
    for j in coe:
        coe2.append(mpf(str(j)))
    
    res = polyroots(coe2, maxsteps=100, extraprec=110)
    GP.append(res)

    dPoly = LkDiff[i]
    weights = []
    for i in res:
        dPolyVal = dPoly.subs(x,i)
        
        ww = mpf('2')/(mpf('1')-i*i)/dPolyVal/dPolyVal
        weights.append(ww)
    
    Weights.append(weights)


temp = "std::vector<std::vector<prec>> GaussPoints1D::xi = {\n"
tempWi = "std::vector<std::vector<prec>> GaussPoints1D::wi = {\n"
for i in GP:
    temp += "{"
    temp += str(i[0])
    for j in range(1,len(i)):
        temp += ", " + str(i[j])
    temp += "},\n"

for i in Weights:
    tempWi += "{"
    tempWi += str(i[0])
    for j in range(1,len(i)):
        tempWi += ", " + str(i[j])
    tempWi += "},\n"

tempWi = tempWi.rstrip(', \n')
tempWi += "\n};"

temp = temp.rstrip(', \n')
temp += "\n};"



headertxt = """
#pragma once

#include "forwarddeclaration.h"
#include "IntegrationPointsBase.h"

#include <vector>


namespace HierAMuS {

 
class GaussPoints1D : public IntegrationPointsBase {
private:
  static std::vector<std::vector<prec>> xi;
  static std::vector<std::vector<prec>> wi;
  indexType maxGP = 0;
  bool init;

public:
  GaussPoints1D();
  ~GaussPoints1D();

  prec getXi(indexType anz, indexType num);
  prec getWeight(indexType anz, indexType num);

  indexType getMaxGP() { return this->xi.size(); };
};

} // namespace HierAMuS
"""

cpptxt = """



#include"GaussPoints1D.h"
#include "control/HandlingStructs.h"


namespace HierAMuS {


GaussPoints1D::GaussPoints1D() {
  this->init = false;


}


GaussPoints1D::~GaussPoints1D() {
}


prec GaussPoints1D::getXi(indexType anz, indexType num) {
  return this->xi[anz - 1][num];
}


prec GaussPoints1D::getWeight(indexType anz, indexType num) {
  return this->wi[anz - 1][num];
}


"""


cpptxt += temp + "\n\n"
cpptxt += tempWi + "\n\n"

cpptxt += """

} // End Namespace

"""


import os, sys

print('sys.argv[0] =', sys.argv[0])             
pathname = os.path.dirname(sys.argv[0])        
print('path =', pathname)
print('full path =', os.path.abspath(pathname)) 

headername = pathname + "/GaussPoints1D.h"
cppname = pathname + "/GaussPoints1D.cpp"

hf=open(headername, 'w')
hf.write(headertxt)
hf.close()

cf=open(cppname, 'w')
cf.write(cpptxt)
cf.close()