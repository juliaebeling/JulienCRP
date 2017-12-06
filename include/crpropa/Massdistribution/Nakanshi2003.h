#ifndef CRPROPA_NAKANSHI2003_H
#define CRPROPA_NAKANSHI2003_H 

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Massdistribution/Density.h"

#include <math>

namespace crpropa {
/**
@class Nakanshi03
@brief Distirbution for HI (arXiv:astro-ph/0304338)
*/

class Nakanshi03: public Density {
private:
	double n0 = 0.94;	//Density at center
	double g1 = 0.6;
	double g2 = 0.24;	//weight for exp-function in Radius

public:
	double getDensity(const Vector3d &position) const {
		
		double z = position.z/kpc;
		double planedensity = getPlanedensity(&position);
		double scaleheight = getScaleheight(&position);
		return planedensity*pow(0.5,pow(z/scaleheight,2))
	}

	double getScaleheight(const Vector3d &position) const {
		
		double R = sqrt(pow(position.x,2)+pow(position.y,2));	//Galactic Radius		
		double scaleheight = 1.06*(116.3 +19.3*R+4.1*pow(R,2)-0.05*pow(R,3));
		return scaleheight;
	}

	double getPlanedensity(const Vector3d &position) const{

		double planedensity = n0*(g1*exp(-R/2.4)+g2*exp(-pow((R-9.5)/4.8),2));
		return planedensity;
	}
};

} //namespace crpopa

#endif 	//CRPROPA_NAKANSHI2003_H 
	



