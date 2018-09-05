#include "crpropa/Massdistribution/Cordes.h"

namespace crpropa {



double Cordes::getHIIDensity(const Vector3d &position) const {
	
	
	double n=0;
	
	double z=position.z;
	double R = sqrt(pow(position.x,2)+pow(position.y,2));	//radius in galactic disk
	
	n += 0.025*exp(-fabs(z/kpc)/1)*exp(-pow(R/20/kpc,2));	//galactocentric component
	n += 0.2*exp(-fabs(z/kpc)/0.15)*exp(-pow((R-4*kpc)/2/kpc,2));	//anular component

	// check if density is NAN
	// return 0 instead
	bool NaN = std::isnan(n);
	if(NaN == true){
		KISS_LOG_WARNING
			<< "\nDensity with 'nan' occured: \n"
			<< "postion = " << position << "\n"
			<< "density-model: Cordes 1991 \n"
			<< "density-type: HII (ionised) \n"
			<< "density is set to 0. \n";
			return 0.;
	}
		return n/ccm;
}

double Cordes::getDensity(const Vector3d &position) const {

	return Cordes::getHIIDensity(position);
}

double Cordes::getNucleonDensity(const Vector3d &position) const {
	
	return getHIIDensity(position);
}

bool Cordes::getisforHI() {
	return isforHI;
}

bool Cordes::getisforHII() {
	return isforHII;
}

bool Cordes::getisforH2() {
	return isforH2;
}

std::string Cordes::getDescription() {
	
	std::stringstream s;
	s << "Density Cordes include HII component";
	return s.str();
}

}//namespace
