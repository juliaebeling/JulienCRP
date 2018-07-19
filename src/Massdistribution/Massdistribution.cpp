#include "crpropa/Massdistribution/Massdistribution.h"


namespace crpropa {

Massdistribution::Massdistribution() {
}

void Massdistribution::add(ref_ptr<Density> dens) { 
			
	bool HI = dens->getisforHI();
	bool HII= dens->getisforHII();
	bool H2 = dens->getisforH2();	
	
	if(HI == true){
		HIDist = dens;
		isforHI = true;
	}
	
	if(HII == true){
		HIIDist = dens;
		isforHII = true;
	}
	
	if(H2 == true){
		H2Dist = dens;
		isforH2 = H2;
	}
	
	return;
	
}


double Massdistribution::getDensity(const Vector3d &position) const{
	
	double n=0.;
	bool noOptionSelected = true;
	
	if(isforHI){
		n = HIDist->getHIDensity(position);
		noOptionSelected = false;
	}
	if(isforHII){
		n += this->HIIDist->getHIIDensity(position);
		noOptionSelected = false;
	}
	if(isforH2){
		n += this->H2Dist->getH2Density(position);
		noOptionSelected = false;
	}
	
	// warning if no Option is load in
	if(noOptionSelected == true){
		KISS_LOG_WARNING
			<< "\ntryed to get density in Massdistribution without loading any Option in or all components are deaktivated.\n"
			<< "Returned density of 0 \n"
			<< "Please load a density in or use a constant density of 0!\n";
			return 0;
	}
	
	return n;
}

double Massdistribution::getHIDensity(const Vector3d &position) const {
	return HIDist->getHIDensity(position);
}

double Massdistribution::getHIIDensity(const Vector3d &position) const {
	return HIIDist->getHIIDensity(position);
}

double Massdistribution::getH2Density(const Vector3d &position) const {
	return H2Dist->getH2Density(position);
}

bool Massdistribution::getisforHI() {
	return isforHI;
}

bool Massdistribution::getisforHII() {
	return isforHII;
}

bool Massdistribution::getisforH2() {
	return isforH2;
}

void Massdistribution::deaktivateHI() {
	isforHI=false;
	return;
}

void Massdistribution::deaktivateHII() {
	isforHII=false;
	return;
}

void Massdistribution::deaktivateH2() {
	isforH2=false;
	return;
}

void MassdistributionSuperposition::addDensity(ref_ptr<Density> dens) {
	DensityList.push_back(dens);
}

double MassdistributionSuperposition::getDensity(const Vector3d &position) const {
	double n = 0.;
	for (int i = 0; i < DensityList.size(); i++)
		n += DensityList[i]->getDensity(position);
	return n;
}

double MassdistributionSuperposition::getHIDensity(const Vector3d &position) const {
	double n = 0.;
	for (int i = 0; i < DensityList.size(); i++)
		n += DensityList[i]->getHIDensity(position);
	return n;
}

double MassdistributionSuperposition::getHIIDensity(const Vector3d &position) const {
	double n = 0.;
	for (int i = 0; i < DensityList.size(); i++)
		n += DensityList[i]->getHIIDensity(position);
	return n;
}

double MassdistributionSuperposition::getH2Density(const Vector3d &position) const {
	double n = 0.;
	for (int i = 0; i < DensityList.size(); i++)
		n += DensityList[i]->getH2Density(position);
	return n;
}

} //namespace crpropa

