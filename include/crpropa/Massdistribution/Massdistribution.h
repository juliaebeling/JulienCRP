#ifndef CRPROPA_MASSDISTRIBUTION_H
#define CRPROPA_MASSDISTRIBUTION_H


#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Massdistribution/Density.h"

#include <string>
#include <vector>

#include "kiss/logger.h"

namespace crpropa {
/** 
 @class Massdistribution
 @brief Density class for superposition of one HI, one HII and one H2 density of different Types. 
No overlaying of more than one density of each type is posible. 
The add function only load activ components in list. 
*/
class Massdistribution: public Density {


ref_ptr<Density> HIDist;
ref_ptr<Density> HIIDist;
ref_ptr<Density> H2Dist;

bool isforHI=false;
bool isforHII=false;
bool isforH2=false;

bool HIisload=false;
bool HIIisload=false;
bool H2isload=false;


public:
	Massdistribution();
	double getDensity(const Vector3d &position) const;
	double getHIDensity(const Vector3d &position) const;
	double getHIIDensity(const Vector3d &position) const;
	double getH2Density(const Vector3d &position) const;
	double getNucleonDensity(const Vector3d &position) const;

	void add(ref_ptr<crpropa::Density> dens);	
	
	bool getisforHI();
	bool getisforHII();
	bool getisforH2();
	
	void deaktivateHI();
	void deaktivateHII();
	void deaktivateH2();

};
/**
 @class MassdistributionSuperposition
 @brief Superposition of density models. 
 the addDensity function adds a new density to the list. 
 The getDensity function cares about acitvated types in loaded densitys. The get(typ)Density doesn't care.
*/

class MassdistributionSuperposition: public Density {

std::vector<ref_ptr<Density>> DensityList ;

public:
	void addDensity(ref_ptr<Density> density);
	double getDensity(const Vector3d &position) const;
	double getHIDensity(const Vector3d &position) const;
	double getHIIDensity(const Vector3d &position) const;
	double getH2Density(const Vector3d &position) const;
};

	
} //namespace crpropa

#endif //CRPROPA_MASSDISTRIBUTION_H


