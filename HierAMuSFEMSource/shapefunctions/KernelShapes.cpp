// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


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
  case 0:{
{
prec t0 = 1;
shape = -2.4494897427831780982*t0;
shapeDeriv = 0;
}
    }break;
  case 1:{
{
prec t0 = 1;
prec t1 = t0*coor;
shape = -3.1622776601683793320*t1;
shapeDeriv = -3.1622776601683793320*t0;
}
    }break;
  case 2:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
shape = 0.93541434669348534640*t0-4.6770717334674267320*t2;
shapeDeriv = -9.3541434669348534640*t1;
}
    }break;
  case 3:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
shape = 3.1819805153394638598*t1-7.4246212024587490062*t3;
shapeDeriv = 3.1819805153394638598*t0-22.273863607376247019*t2;
}
    }break;
  case 4:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
shape = -0.58630196997792869432*t0+8.2082275796910017205*t2-12.312341369536502581*t4;
shapeDeriv = 16.416455159382003441*t1-49.249365478146010323*t3;
}
    }break;
  case 5:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
shape = -3.1868871959954905188*t1+19.121323175972943113*t3-21.033455493570237424*t5;
shapeDeriv = -3.1868871959954905188*t0+57.363969527918829338*t2-105.16727746785118712*t4;
}
    }break;
  case 6:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
shape = 0.42790824805091102614*t0-11.553522697374597706*t2+42.362916557040191588*t4-36.714527682768166043*t6;
shapeDeriv = -23.107045394749195411*t1+169.45166622816076635*t3-220.28716609660899626*t5;
}
    }break;
  case 7:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
shape = 3.1888018174935236950*t1-35.076819992428760645*t3+91.199731980314777677*t5-65.142665700224841198*t7;
shapeDeriv = 3.1888018174935236950*t0-105.23045997728628194*t2+455.99865990157388839*t4-455.99865990157388839*t6;
}
    }break;
  case 8:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
shape = -0.33711639078736589962*t0+14.833121194644099583*t2-96.415287765186647292*t4+192.83057553037329458*t6-117.07570657201235743*t8;
shapeDeriv = 29.666242389288199167*t1-385.66115106074658917*t3+1156.9834531822397675*t5-936.60565257609885941*t7;
}
    }break;
  case 9:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
shape = -3.1897395624976187074*t1+55.288819083292057595*t3-248.79968587481425918*t5+402.81853903541356248*t7-212.59867337980160242*t9;
shapeDeriv = -3.1897395624976187074*t0+165.86645724987617279*t2-1243.9984293740712959*t4+2819.7297732478949374*t6-1913.3880604182144218*t8;
}
    }break;
  case 10:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
shape = 0.27818150321412232602*t0-18.081797708917951191*t2+180.81797708917951191*t4-614.78112210321034050*t6+834.34580856864260496*t8-389.36137733203321565*t10;
shapeDeriv = -36.163595417835902382*t1+723.27190835671804764*t3-3688.6867326192620430*t5+6674.7664685491408397*t7-3893.6137733203321565*t9;
}
    }break;
  case 11:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
prec t11 = t10*coor;
shape = 3.1902669229314937136*t1-79.756673073287342840*t3+542.34537689835393131*t5-1472.0803087241035278*t7+1717.4270268447874492*t9-718.19675668054747874*t11;
shapeDeriv = 3.1902669229314937136*t0-239.27001921986202852*t2+2711.7268844917696566*t4-10304.562161068724695*t6+15456.843241603087042*t8-7900.1643234860222661*t10;
}
    }break;
  case 12:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
prec t11 = t10*coor;
prec t12 = t11*coor;
shape = -0.23681590286673303879*t0+21.313431258005973491*t2-301.94027615508462446*t4+1529.8307325190954306*t6-3442.1191481679647188*t8+3518.6106847939194903*t10-1332.8070775734543524*t12;
shapeDeriv = 42.626862516011946982*t1-1207.7611046203384978*t3+9178.9843951145725835*t5-27536.953185343717751*t7+35186.106847939194903*t9-15993.684930881452229*t11;
}
    }break;
  case 13:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
prec t11 = t10*coor;
prec t12 = t11*coor;
prec t13 = t12*coor;
shape = -3.1905924437652506392*t1+108.48014308801852173*t3-1030.5613593361759565*t5+4122.2454373447038259*t7-7900.9704215773489996*t9+7182.7003832521354541*t11-2486.3193634334315034*t13;
shapeDeriv = -3.1905924437652506392*t0+325.44042926405556520*t2-5152.8067966808797823*t4+28855.718061412926781*t6-71108.733794196140996*t8+79009.704215773489996*t10-32322.151724634609544*t12;
}
    }break;
  case 14:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
prec t11 = t10*coor;
prec t12 = t11*coor;
prec t13 = t12*coor;
prec t14 = t13*coor;
shape = 0.20617366808783367477*t0-24.534666502452207298*t2+466.15866354659193866*t4-3263.1106448261435706*t6+10721.649261571614589*t8-17869.415435952690982*t10+14620.430811234019894*t12-4659.2581706130393070*t14;
shapeDeriv = -49.069333004904414595*t1+1864.6346541863677546*t3-19578.663868956861424*t5+85773.194092572916713*t7-178694.15435952690982*t9+175445.16973480823873*t11-65229.614388582550297*t13;
}
    }break;
  case 15:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
prec t11 = t10*coor;
prec t12 = t11*coor;
prec t13 = t12*coor;
prec t14 = t13*coor;
prec t15 = t14*coor;
shape = 3.1908073201802004956*t1-141.45912452798888864*t3+1782.3849690526599968*t5-9760.6795924312333160*t7+27112.998867864536989*t9-39930.052878127772656*t11+29691.577781171933514*t13-8766.0848687269517993*t15;
shapeDeriv = 3.1908073201802004956*t0-424.37737358396666591*t2+8911.9248452632999841*t4-68324.757147018633212*t6+244016.98981078083290*t8-439230.58165940549922*t10+385990.51115523513568*t12-131491.27303090427699*t14;
}
    }break;
  case 16:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
prec t11 = t10*coor;
prec t12 = t11*coor;
prec t13 = t12*coor;
prec t14 = t13*coor;
prec t15 = t14*coor;
prec t16 = t15*coor;
shape = -0.18255978972530029810*t0+27.749088038245645311*t2-679.85265693701831011*t4+6254.6444438205684530*t6-27922.519838484680594*t8+67014.047612363233425*t10-88336.699125387898606*t12+60185.443360154392457*t14-16550.996924042457926*t16;
shapeDeriv = 55.498176076491290621*t1-2719.4106277480732405*t3+37527.866662923410718*t5-223380.15870787744475*t7+670140.47612363233425*t9-1060040.3895046547833*t11+842596.20704216149440*t13-264815.95078467932681*t15;
}
    }break;
  case 17:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
prec t11 = t10*coor;
prec t12 = t11*coor;
prec t13 = t12*coor;
prec t14 = t13*coor;
prec t15 = t14*coor;
prec t16 = t15*coor;
prec t17 = t16*coor;
shape = -3.1909565313996316043*t1+178.69356575837936984*t3-2876.9664087099078545*t5+20549.760062213627532*t7-77061.600233301103244*t9+162529.92049205323593*t11-193785.67443283270438*t13+121808.13821492341418*t15-31347.682628840584532*t17;
shapeDeriv = -3.1909565313996316043*t0+536.08069727513810953*t2-14384.832043549539272*t4+143848.32043549539272*t6-693554.40209970992920*t8+1787829.1254125855953*t10-2519213.7676268251570*t12+1827122.0732238512127*t14-532910.60469028993705*t16;
}
    }break;

default:
    throw std::runtime_error("LegendreShapes::getShape: functionNumber out of range");
}

}

auto KernelShapes::getShape(prec coor, indexType functionNumber) -> KernelShapesValues {

  KernelShapesValues shape;
  switch(functionNumber) {
  case 0:{
{
prec t0 = 1;
shape.shapeValue = -2.4494897427831780982*t0;
shape.shapeDerivative = 0;
}
    }break;
  case 1:{
{
prec t0 = 1;
prec t1 = t0*coor;
shape.shapeValue = -3.1622776601683793320*t1;
shape.shapeDerivative = -3.1622776601683793320*t0;
}
    }break;
  case 2:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
shape.shapeValue = 0.93541434669348534640*t0-4.6770717334674267320*t2;
shape.shapeDerivative = -9.3541434669348534640*t1;
}
    }break;
  case 3:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
shape.shapeValue = 3.1819805153394638598*t1-7.4246212024587490062*t3;
shape.shapeDerivative = 3.1819805153394638598*t0-22.273863607376247019*t2;
}
    }break;
  case 4:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
shape.shapeValue = -0.58630196997792869432*t0+8.2082275796910017205*t2-12.312341369536502581*t4;
shape.shapeDerivative = 16.416455159382003441*t1-49.249365478146010323*t3;
}
    }break;
  case 5:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
shape.shapeValue = -3.1868871959954905188*t1+19.121323175972943113*t3-21.033455493570237424*t5;
shape.shapeDerivative = -3.1868871959954905188*t0+57.363969527918829338*t2-105.16727746785118712*t4;
}
    }break;
  case 6:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
shape.shapeValue = 0.42790824805091102614*t0-11.553522697374597706*t2+42.362916557040191588*t4-36.714527682768166043*t6;
shape.shapeDerivative = -23.107045394749195411*t1+169.45166622816076635*t3-220.28716609660899626*t5;
}
    }break;
  case 7:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
shape.shapeValue = 3.1888018174935236950*t1-35.076819992428760645*t3+91.199731980314777677*t5-65.142665700224841198*t7;
shape.shapeDerivative = 3.1888018174935236950*t0-105.23045997728628194*t2+455.99865990157388839*t4-455.99865990157388839*t6;
}
    }break;
  case 8:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
shape.shapeValue = -0.33711639078736589962*t0+14.833121194644099583*t2-96.415287765186647292*t4+192.83057553037329458*t6-117.07570657201235743*t8;
shape.shapeDerivative = 29.666242389288199167*t1-385.66115106074658917*t3+1156.9834531822397675*t5-936.60565257609885941*t7;
}
    }break;
  case 9:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
shape.shapeValue = -3.1897395624976187074*t1+55.288819083292057595*t3-248.79968587481425918*t5+402.81853903541356248*t7-212.59867337980160242*t9;
shape.shapeDerivative = -3.1897395624976187074*t0+165.86645724987617279*t2-1243.9984293740712959*t4+2819.7297732478949374*t6-1913.3880604182144218*t8;
}
    }break;
  case 10:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
shape.shapeValue = 0.27818150321412232602*t0-18.081797708917951191*t2+180.81797708917951191*t4-614.78112210321034050*t6+834.34580856864260496*t8-389.36137733203321565*t10;
shape.shapeDerivative = -36.163595417835902382*t1+723.27190835671804764*t3-3688.6867326192620430*t5+6674.7664685491408397*t7-3893.6137733203321565*t9;
}
    }break;
  case 11:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
prec t11 = t10*coor;
shape.shapeValue = 3.1902669229314937136*t1-79.756673073287342840*t3+542.34537689835393131*t5-1472.0803087241035278*t7+1717.4270268447874492*t9-718.19675668054747874*t11;
shape.shapeDerivative = 3.1902669229314937136*t0-239.27001921986202852*t2+2711.7268844917696566*t4-10304.562161068724695*t6+15456.843241603087042*t8-7900.1643234860222661*t10;
}
    }break;
  case 12:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
prec t11 = t10*coor;
prec t12 = t11*coor;
shape.shapeValue = -0.23681590286673303879*t0+21.313431258005973491*t2-301.94027615508462446*t4+1529.8307325190954306*t6-3442.1191481679647188*t8+3518.6106847939194903*t10-1332.8070775734543524*t12;
shape.shapeDerivative = 42.626862516011946982*t1-1207.7611046203384978*t3+9178.9843951145725835*t5-27536.953185343717751*t7+35186.106847939194903*t9-15993.684930881452229*t11;
}
    }break;
  case 13:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
prec t11 = t10*coor;
prec t12 = t11*coor;
prec t13 = t12*coor;
shape.shapeValue = -3.1905924437652506392*t1+108.48014308801852173*t3-1030.5613593361759565*t5+4122.2454373447038259*t7-7900.9704215773489996*t9+7182.7003832521354541*t11-2486.3193634334315034*t13;
shape.shapeDerivative = -3.1905924437652506392*t0+325.44042926405556520*t2-5152.8067966808797823*t4+28855.718061412926781*t6-71108.733794196140996*t8+79009.704215773489996*t10-32322.151724634609544*t12;
}
    }break;
  case 14:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
prec t11 = t10*coor;
prec t12 = t11*coor;
prec t13 = t12*coor;
prec t14 = t13*coor;
shape.shapeValue = 0.20617366808783367477*t0-24.534666502452207298*t2+466.15866354659193866*t4-3263.1106448261435706*t6+10721.649261571614589*t8-17869.415435952690982*t10+14620.430811234019894*t12-4659.2581706130393070*t14;
shape.shapeDerivative = -49.069333004904414595*t1+1864.6346541863677546*t3-19578.663868956861424*t5+85773.194092572916713*t7-178694.15435952690982*t9+175445.16973480823873*t11-65229.614388582550297*t13;
}
    }break;
  case 15:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
prec t11 = t10*coor;
prec t12 = t11*coor;
prec t13 = t12*coor;
prec t14 = t13*coor;
prec t15 = t14*coor;
shape.shapeValue = 3.1908073201802004956*t1-141.45912452798888864*t3+1782.3849690526599968*t5-9760.6795924312333160*t7+27112.998867864536989*t9-39930.052878127772656*t11+29691.577781171933514*t13-8766.0848687269517993*t15;
shape.shapeDerivative = 3.1908073201802004956*t0-424.37737358396666591*t2+8911.9248452632999841*t4-68324.757147018633212*t6+244016.98981078083290*t8-439230.58165940549922*t10+385990.51115523513568*t12-131491.27303090427699*t14;
}
    }break;
  case 16:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
prec t11 = t10*coor;
prec t12 = t11*coor;
prec t13 = t12*coor;
prec t14 = t13*coor;
prec t15 = t14*coor;
prec t16 = t15*coor;
shape.shapeValue = -0.18255978972530029810*t0+27.749088038245645311*t2-679.85265693701831011*t4+6254.6444438205684530*t6-27922.519838484680594*t8+67014.047612363233425*t10-88336.699125387898606*t12+60185.443360154392457*t14-16550.996924042457926*t16;
shape.shapeDerivative = 55.498176076491290621*t1-2719.4106277480732405*t3+37527.866662923410718*t5-223380.15870787744475*t7+670140.47612363233425*t9-1060040.3895046547833*t11+842596.20704216149440*t13-264815.95078467932681*t15;
}
    }break;
  case 17:{
{
prec t0 = 1;
prec t1 = t0*coor;
prec t2 = t1*coor;
prec t3 = t2*coor;
prec t4 = t3*coor;
prec t5 = t4*coor;
prec t6 = t5*coor;
prec t7 = t6*coor;
prec t8 = t7*coor;
prec t9 = t8*coor;
prec t10 = t9*coor;
prec t11 = t10*coor;
prec t12 = t11*coor;
prec t13 = t12*coor;
prec t14 = t13*coor;
prec t15 = t14*coor;
prec t16 = t15*coor;
prec t17 = t16*coor;
shape.shapeValue = -3.1909565313996316043*t1+178.69356575837936984*t3-2876.9664087099078545*t5+20549.760062213627532*t7-77061.600233301103244*t9+162529.92049205323593*t11-193785.67443283270438*t13+121808.13821492341418*t15-31347.682628840584532*t17;
shape.shapeDerivative = -3.1909565313996316043*t0+536.08069727513810953*t2-14384.832043549539272*t4+143848.32043549539272*t6-693554.40209970992920*t8+1787829.1254125855953*t10-2519213.7676268251570*t12+1827122.0732238512127*t14-532910.60469028993705*t16;
}
    }break;

default:
    throw std::runtime_error("LegendreShapes::getShape: functionNumber out of range");
}

  return shape;
}

}
