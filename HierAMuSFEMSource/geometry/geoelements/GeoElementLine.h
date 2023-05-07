// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#ifndef GEOELEMENTLINE_H_
#define GEOELEMENTLINE_H_

#define GEOLTEMP 
#define GEONAME GeoElementLine::


#include <GeoElement.h>
#include <Edge.h>
#include <iostream>

namespace HierAMuS {

GEOLTEMP
class GeoElementLine: public GeoElement {
public:
	GeoElementLine();
	~GeoElementLine();
	void setVertex(std::vector<Vertex*> &in);
private:
	Edge *edgePtr;
	Vertex *vert1, *vert2;
};

GEOLTEMP
GEONAME GeoElementLine(){
	this->edgePtr = 0;
	this->vert1 = 0;
	this->vert2 = 0;
}

GEOLTEMP
GEONAME ~GeoElementLine(){

}

GEOLTEMP
void GEONAME setVertex(std::vector<Vertex*> &in){
	if(in.size()!=2){

	}else{
		this->vert1 = in[0];
		this->vert2 = in[1];
		std::cout << "Richtig" << std::endl;
	}
}

} /* namespace HierAMuS */

#undef GEOTEMP
#undef GEONAME

#endif /* GEOELEMENTLINE_H_ */
