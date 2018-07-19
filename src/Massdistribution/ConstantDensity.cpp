#include "crpropa/Massdistribution/ConstantDensity.h"

namespace crpropa{

constantDensity::constantDensity(double HI, double HII, double H2) {
	if(HI!=0)
		setHI(true, HI);
	if(HII!=0)
		setHII(true, HII);
	if(H2!=0)
		setH2(true, H2);
}

constantDensity::constantDensity(double densitynumber) {
	if(densitynumber == 0)
	{
		setHI(true,0);
	}
	else
	{
		setHI(true,densitynumber);
		setHII(true,densitynumber);
		setH2(true,densitynumber);
	}
}

double constantDensity::getDensity(const Vector3d &position) const {
	double n = 0;
				
	if(isforHI) 
		n += HIdensitynumber;
	if(isforHII)
		n += HIIdensitynumber;
	if(isforH2)
		n += H2densitynumber;
	
	//check if any density is activ and give warning if not
	bool anyDensityActive = isforHI||isforHII||isforH2;

	if(anyDensityActive == false){
		KISS_LOG_WARNING
			<< "\n tryed to get density although all density-types are deaktivated \n"
			<< "density-module: constantDensity\n"
			<< "returned 0 density\n"
			<< "please use constant Density with 0 \n";
	}

	return n;
}

double constantDensity::getHIDensity(const Vector3d &position) const {
		
	return HIdensitynumber;
}
	
double constantDensity::getHIIDensity(const Vector3d &position) const{
		
	return HIIdensitynumber;
}

double constantDensity::getH2Density(const Vector3d &position) const{

	return H2densitynumber;
}

bool constantDensity::getisforHI() {

	return isforHI;
}

bool constantDensity::getisforHII() {

	return isforHII;
}

bool constantDensity::getisforH2() {

	return isforH2;
}

void constantDensity::setHI(bool activate, double densitynumber) {	
	
	isforHI = activate;
	HIdensitynumber = densitynumber;
}

void constantDensity::setHI(bool activate) {
	setHI(activate, HIdensitynumber);	 
}

void constantDensity::setHI(double densitynumber) {
	setHI(isforHI, densitynumber);	
}


void constantDensity::setHII(bool activate, double densitynumber) {

	isforHII = activate;
	HIIdensitynumber = densitynumber;
}

void constantDensity::setHII(bool activate) {
	setHII(activate, HIIdensitynumber);	 
}

void constantDensity::setHII(double densitynumber) {
	setHII(isforHII, densitynumber);	
}

void constantDensity::setH2(bool activate, double densitynumber) {

	isforH2 = activate;
	H2densitynumber = densitynumber;
}
void constantDensity::setH2(bool activate) {
	setH2(activate, H2densitynumber);	 
}

void constantDensity::setH2(double densitynumber) {
	setH2(isforH2, densitynumber);	
}

}//namespace 
