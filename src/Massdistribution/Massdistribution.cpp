#include "crpropa/Massdistribution/Massdistribution.h"


namespace crpropa {

void Massdistribution::add(ref_ptr<Density> dens) { 
			
	bool HI = dens->getisforHI();	// check which density tpye is activated in loading density. Just use this part!
	bool HII= dens->getisforHII();
	bool H2 = dens->getisforH2();	
	
	bool nothingToLoad = !(HI || HII || H2);
	
	if(nothingToLoad == true){
		KISS_LOG_WARNING
		<<"\n tryed to add density to Massdistribution wether no density type is activated. \n  nothing is load!\n";
		return ;
	}
	
	if(HI == true){
		HIDist = dens;
		isforHI = true;
		HIisload = true;
	}
	
	if(HII == true){
		HIIDist = dens;
		isforHII = true;
		HIIisload = true;
	}
	
	if(H2 == true){
		H2Dist = dens;
		isforH2 = H2;
		H2isload = true;
	}
	
	return;
	
}


double Massdistribution::getDensity(const Vector3d &position) const{
	
	double n=0.;

	// warning if nothing is load in 	
	bool nothingLoadIn = !(HIisload || HIIisload || H2isload);
	if(nothingLoadIn)
	{	
		KISS_LOG_WARNING
		<<"\n tryed to get density in Massdistribution allthough no densitytype is load in \n"
		<<" density return is 0 \n";
		return 0;
	} 
	
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
	
	// warning if no option is activ in
	if(noOptionSelected == true){
		KISS_LOG_WARNING
			<< "\ntryed to get density in Massdistribution where all components are deaktivated.\n"
			<< "Returned density of 0 \n"
			<< "Please activate a density in or use a constant density of 0!\n";
			return 0;
	}
	
	return n;
}

double Massdistribution::getNucleonDensity(const Vector3d &position) const{
	
	double n=0.;
	bool nothingLoadIn = !(HIisload || HIIisload || H2isload);
	
	if(nothingLoadIn)
	{	
		KISS_LOG_WARNING
		<<"\n tryed to get NucleonDensity in Massdistribution allthough no densitytype is load in \n"
		<<" density return is 0 \n";
		return 0;
	} 

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
		n += 2*this->H2Dist->getH2Density(position);
		noOptionSelected = false;
	}
	
	// warning if no option is activ
	if(noOptionSelected == true){
		KISS_LOG_WARNING
			<< "\ntryed to get nucleon-density in Massdistribution without loading any Option in or all components are deaktivated.\n"
			<< "Returned density of 0 \n"
			<< "Please load a density in or use a constant density of 0!\n";
			return 0;
	}
	
	return n;
}

double Massdistribution::getHIDensity(const Vector3d &position) const {
	
	if(HIisload)
		return HIDist->getHIDensity(position);
	
	// warning if no HI is load
	KISS_LOG_WARNING 
	<< "\n tryed to get HIDensity in Massdistribution, allthough no HI option is load in. \n"
	<< "Please load HI Option or use constant Density of 0 \n"
	<< "return a density of 0.\n";
	return 0.;
}

double Massdistribution::getHIIDensity(const Vector3d &position) const {
	
	if(HIIisload)
		return HIIDist->getHIIDensity(position);
	
	// warning if no HII is load
	KISS_LOG_WARNING 
	<< "\n tryed to get HIIDensity in Massdistribution, allthough no HII option is load in. \n"
	<< "Please load HII Option or use constant Density of 0 \n"
	<< "return a density of 0.\n";
	return 0.;
}

double Massdistribution::getH2Density(const Vector3d &position) const {
		
	if(H2isload)
		return H2Dist->getH2Density(position);
	
	// warning if no H2 is load
	KISS_LOG_WARNING 
	<< "\n tryed to get H2Density in Massdistribution, allthough no H2 option is load in. \n"
	<< "Please load H2 Option or use constant Density of 0 \n"
	<< "return a density of 0.\n";
	return 0.;
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

double MassdistributionSuperposition::getNucleonDensity(const Vector3d &position) const {
	double n = 0.;
	for (int i = 0; i < DensityList.size(); i++)
		n += DensityList[i]->getNucleonDensity(position);
	return n;
}

} //namespace crpropa

