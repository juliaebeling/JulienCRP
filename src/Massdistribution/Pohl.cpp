#include "crpropa/Massdistribution/Pohl.h"
#include "crpropa/Module.h"



#include <fstream>

namespace crpropa {
Pohl::Pohl() {
	loadGridHI();
	loadGridH2();
}

void Pohl::loadGridHI() {

}

void Pohl::loadGridH2() {


}

double Pohl::getH2Density(const Vector3d &position) const {
	
	double n = 0; //density in ccm
	
	Vector3d pos = position;
	
	// set galactocentric coordinate system with the Sun at (-8.5,0.,0.) instead of (8.5, 0, 0) to be consistand with JF12 implementation
	pos.x = -position.x;
	pos.y = -position.y;
	
	
	if(fabs(pos.x)>15*kpc)	//boundary of Grid not repeat
	{
		return 0;
	}
	if(fabs(pos.y)>15*kpc)
	{
		return 0;
	}
	if(fabs(pos.z)>500*pc)
	{
		return 0;
	}
	
	n= H2density.interpolate(pos);
		
	return n/ccm;
}


double Pohl::getHIDensity(const Vector3d &position) const {
	
	double n = 0;	// density in ccm
	Vector3d pos = position;
	
	// set galactocentric coordinate system with the Sun at (-8.5,0.,0.) instead of (8.5, 0, 0) to be consistand with JF12 implementation
	pos.x = -position.x;
	pos.y = -position.y;
	
	if(fabs(pos.x)>20*kpc)	//boundary of Grid not repeat 
	{
		return 0;
	}
	if(fabs(pos.y)>20*kpc)
	{
		return 0;
	}
	if(fabs(pos.z)>1500*pc)
	{
		return 0;
	}
	
	n= HIdensity.interpolate(pos);
	
	
	return n/ccm;
	
}
double Pohl::getDensity(const Vector3d &position) const {
	double n=0;
	if(isforHI)
		n+=Pohl::getHIDensity(position);
	if(isforH2)
		n+=Pohl::getH2Density(position);
		
	//check if any density is activ and give warning if not
	bool anyDensityActive = isforHI||isforH2;

	if(anyDensityActive == false){
		KISS_LOG_WARNING
			<< "\n called getDensity on deactivated Pohl \n"
			<< "returned 0 density\n"
			<< "please activate \n";
	}	
	
	return n;
}


double Pohl::getNucleonDensity(const Vector3d &position) const {
	double n=0;
	if(isforHI)
		n+= Pohl::getHIDensity(position);
	if(isforH2)
		n+= 2*Pohl::getH2Density(position);
		
	//check if any density is activ and give warning if not
	bool anyDensityActive = isforHI||isforH2;

	if(anyDensityActive == false){
		KISS_LOG_WARNING
			<< "\n called getNucleonDensity on deactivated Ferriere \n"
			<< "returned 0 density\n"
			<< "please activate \n"; 
	}	
	
	return n;
}


bool Pohl::getisforHI() {
	return isforHI;
}

bool Pohl::getisforHII() {
	return isforHII;
}

bool Pohl::getisforH2() {
	return isforH2;
}

void Pohl::setisforHI(bool HI) {
	isforHI=HI;
}

void Pohl::setisforH2(bool H2) {
	isforH2=H2;
}

std::string Pohl::getDescription() {
	
	std::stringstream s;
	s << "Density modell Pohl 2008: ";
	s<< "HI component is ";
	if(isforHI==false)
		s<< "not ";
	s<< "activ. HII component is ";
	if(isforH2==false)
		s<<"not "; 
	s<<"activ. Pohl has no HII component.";
	
	return s.str();
}

} //namespace
