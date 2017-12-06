#ifndef CRPROPA_MASSDISTRIBUTION_H
#define CRPROPA_MASSDISTRIBUTION_H


#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Massdistribution/Density.h"
#include "crpropa/Massdistribution/Nakanshi2003.h"
#include "crpropa/Massdistribution/Nakanshi2006.h"



namespace crpropa {

class Masssdistribution: public Density {
	std::vector<ref_ptr<Density> >Distributions;
	int HIFlag;		// 0: none	1: Nakanshi03	2:Ferrie07
	int HIIFlag;		// 0: none	1: Cordes91	2:Ferrie07	
	int H2Flag;		// 0: none	1: Nakanshi06	2:Ferrie07	3:Pohl08
	
public:
	Massdistribution();
	Massdistribution(int HI, int HII, int H2);
	double getDensity(const vector3d &position) const;
	void setHI(int HI);
	void setHII(int HII);
	void setH2(int H2);
	void setAll(int HI, int HII, int H2);
	
} //namespace crpropa

#endif //CRPROPA_MASSDISTRIBUTION_H


