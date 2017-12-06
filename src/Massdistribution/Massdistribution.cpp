#include "crpropa/Massdistribution/Massdistribution.h"


namespace crpropa {
	
Massdistribution::Massdistribution() {
	setAll(0,0,0);
}

Massdistribution::Massdistribution(int HI, int HII, int H2){
	setAll(HI,HII,H2);
}

double Massdistribution::getDensity(const vector3d &position) const{
	
	double n
	for(int i=0; i < Distributions.size(); i++)
		n += Distritbutions[i]->getDensity(position);
	return b;
}

void Massdistribution::setHI(int HI){

	this->HI = HI;
	switch(HI){
	case 0:
		Distributions[0] = new Destiny();
	break;
	case 1: Distributions[0] = new Nakanshi03(); 
	break;
	default: 
		throw std::runtime_error("Massdistribution: unexpected Distribution for HI");
	}
	return;
}

void Massdistribution::setHII(int HII){

	this->HII = HII;
	switch(HII){
	case 0:
		Distribution[1] = new Density();
	break;
	case 1:
		Distribution[1] = new Nakanshi06();
	default:
		throw std::runtime_error("Massdistribution: unexpected Distribution for HII");
	}
	return;
}

void Massdistribution::setH2(int H2){

	this->H2 = H2;
	switch(H2){
	case 0:
		Distribution[2] = new Density();
	break;
	default: 
		throw std::runtime_error("Massdistribution: unexpected Distribution for H2");
	}
	return;
}

void Massdistribution::setAll(int HI, int HII, int H2){
	
	setHI(HI);
	setHII(HII);
	setH2(H2);
	return;
}

} //namespace crpropa

