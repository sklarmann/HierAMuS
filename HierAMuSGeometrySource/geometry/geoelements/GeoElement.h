// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#ifndef GEOELEMENT_H_
#define GEOELEMENT_H_

#define GEOTEMP 

#include <vector>
#include <Vertex.h>

namespace HierAMuS {

GEOTEMP
class GeoElement {
public:
	GeoElement();
	virtual ~GeoElement();
	virtual void setId(indexType &iid) {id=iid;};
	virtual short getNumberOfVertices();
	virtual short getNumberOfEdges();
	virtual short getNumberOfFaces();
	virtual void setVertex(std::vector<Vertex*> &in);
protected:
	indexType id;
};

GEOTEMP
GeoElement::GeoElement(){

}

GEOTEMP
GeoElement::~GeoElement(){

}

GEOTEMP
short GeoElement::getNumberOfVertices(){

	return 0;
}

GEOTEMP
short GeoElement::getNumberOfEdges(){

	return 0;
}


GEOTEMP
short GeoElement::getNumberOfFaces(){
	return 0;
}

GEOTEMP
void GeoElement::setVertex(std::vector<Vertex*> &in){

}
} /* namespace HierAMuS */


#undef GEOTEMP
#endif /* GEOELEMENT_H_ */
