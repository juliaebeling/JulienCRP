#ifndef CRPROPA_NAKANSHI2006_H
#define CRPROPA_NAKANSHI2006_H 

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Massdistribution/Density.h"

#include <math>

namespace crpropa {
/**
@class Nakanshi06
@brief Distirbution for H2 (arxiv: astro-ph/0610769)
*/

class Nakanshi06: public Density {
private:
	double g1 = 11.2;
	double g2 = 0.83;	//weight for radial exp-fkt
	double n0 = 0.94;	//density at center
	double z0 = 42.78;	//offset of schaleheight
	double g3 = 10.8;	//weight of schaleheight exp-fkt
	double zH = 1.06e-3;	//standart scaleheight
public:
	
	Nakanshi06() {
	}

	double getPlaneDensity(const vector3d &position) const {

	double R = sqrt(pow(position.x,2)+pow(position.y,2))/kpc;
	double Planddensity = g1*exp(-pow(R,2)/0.85) +g2*exp(-pow((r-4)/3.2),2));
	return n0*Planedensity;
	}
	
	double getScaleheight(const vector3d &postion) const {

	double R = sqrt(pow(postion.x,2)+ pow(position.y,2))/kpc;
	double Scaleheight = zH*(g1*exp(0.28*R)+z0);
	return Scaleheight;
	}

	double getDensity(const vector3d &postion) const {

	double z = positon.z/kpc
	double plane = getPlaneDensity(&position);
	double scaleheight = getScaleheight(&position);
	return 2*plane*pow(0.5,pow(position.z,2));	//Factor 2 for both nuclei
	}
};

} //namespace

#endif   // CRPROPA_NAKANSHI2006_H
	

