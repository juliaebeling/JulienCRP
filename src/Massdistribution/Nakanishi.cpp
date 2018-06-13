#include "crpropa/Massdistribution/Nakanishi.h"

namespace crpropa{ 

Nakanishi::Nakanishi(){
	isforHI = true;
	isforHII = false;
	isforH2 = true;		
}



double Nakanishi::getHIScaleheight(const Vector3d &position) {
	
	double x = position.x/kpc;
	double y = position.y/kpc;
	double R = sqrt(pow(x,2)+pow(y,2));	//Galactic Radius		
	double scaleheight = 1.06*(116.3 +19.3*R+4.1*pow(R,2)-0.05*pow(R,3));
	return scaleheight*pc;
	}

double Nakanishi::getHIPlanedensity(const Vector3d &position) {
	
	double x = position.x/kpc;
	double y = position.y/kpc;
	double R = sqrt(pow(x,2)+pow(y,2));
	double planedensity = 0.94*(0.6*exp(-R/2.4)+0.24*exp(-pow((R-9.5)/4.8,2)));
	return planedensity;
	}


double Nakanishi::getH2Scaleheight(const Vector3d &position)  {

	double x = position.x/kpc;
	double y = position.y/kpc;
	double R = sqrt(pow(x,2)+ pow(y,2));
	double Scaleheight = 1.06e-3*( 11.2*exp(0.28*R)+42.78);
	return Scaleheight*kpc;
}

double Nakanishi::getH2Planedensity(const Vector3d &position) {

	double x = position.x/kpc;
	double y = position.y/kpc;
	double R = sqrt(pow(x,2)+pow(y,2));
	double Planedensity = 11.2*exp(-pow(R,2)/0.874) +0.83*exp(-pow((R-4)/3.2,2));
	return 0.94*Planedensity;
}

double Nakanishi::getHIDensity(const Vector3d &position) {
	
	double z = position.z/kpc;
	double planedensity = getHIPlanedensity(position);
	double scaleheight = getHIScaleheight(position)/kpc;
	return planedensity*pow(0.5,pow(z/scaleheight,2));
}

double Nakanishi::getH2Density(const Vector3d &position) {
	
	double z = position.z/kpc;
	double plane = getH2Planedensity(position);
	double scaleheight = getH2Scaleheight(position)/kpc;
	return plane*pow(0.5,pow(z/scaleheight,2));
}
	

double Nakanishi::getDensity(const Vector3d &position) {
	double n = 0;
	if(isforHI)
		n += getHIDensity(position);
	if(isforH2)
		n += getH2Density(position);
	return n;
}

bool Nakanishi::getisforHI() {
	return isforHI;
}

bool Nakanishi::getisforHII() {
	return isforHII;
}
bool Nakanishi::getisforH2() {
	return isforH2;
}

void Nakanishi::setisforHI(bool HI) {
	isforHI = HI;
	if(isforHI)
		return;
	if(isforH2)
		return;
	isforH2=true;
}

void Nakanishi::setisforH2(bool H2) {
	isforH2 = H2;
	if(isforHI)
		return;
	if(isforH2)
		return;
	isforH2=true;
}

} //namespace crpropa









