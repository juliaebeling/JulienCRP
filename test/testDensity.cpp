#include "crpropa/Massdistribution/Massdistribution.h"
#include "crpropa/Massdistribution/Cordes.h"
#include "crpropa/Massdistribution/Ferriere07.h"
#include "crpropa/Massdistribution/Nakanishi.h"
#include "crpropa/Massdistribution/Pohl08.h"
#include "crpropa/Massdistribution/ConstantDensity.h"
#include "crpropa/Units.h"

#include <stdexcept>
#include <cmath>

namespace crpropa {

TEST(testConstantDensity, SimpleTest) {
	//test constant Density in all types and in total density (output)
	constantDensity n(2/ccm,3/ccm, 2/ccm);
	Vector3d p(1*pc,2*pc,1*kpc); 	// random position for testing density
	double HI = n.getHIDensity(p);
	double HII = n.getHIIDensity(p);
	double H2 = n.getH2Density(p);
	double Hges = n.getDensity(p);
	EXPECT_DOUBLE_EQ(HI, 2e6);	// density output in m^-3
	EXPECT_DOUBLE_EQ(HII, 3e6);
	EXPECT_DOUBLE_EQ(H2, 2e6);
	EXPECT_DOUBLE_EQ(Hges, 7e6);	// total density 2+3+2 = 7 (/ccm)
	
	//test set/get function for used type
	
	bool useHI = n.getisforHI();
	bool useHII= n.getisforHII();
	bool useH2 = n.getisforH2();
	EXPECT_BOOL_EQ(useHI, true); //all types are activated
	EXPECT_BOOL_EQ(useHII,true);
	EXPECT_BOOL_EQ(useH2, true);
	
	//set density number to 500
	n.setHI(500.);
	n.setHII(500.);
	n.setH2(500.);
	
	HI = n.getHIDensity(p);		
	HII = n.getHIIDensity(p);
	H2 = n.getH2Density(p);
	
	//check if output is changed too 500 in types and is 0 (all deaktivated) for Hges
	EXPECT_DOUBLE_EQ(HI, 500);	
	EXPECT_DOUBLE_EQ(HII, 500);
	EXPECT_DOUBLE_EQ(H2, 500);
	
	//set all type-use to false
	n.setHI(false);
	n.setHII(false);
	n.setH2(false);
	
	//check if all isfor are set to false
	useHI = n.getisforHI();
	useHII= n.getisforHII();
	useH2 = n.getisforH2();
	EXPECT_BOOL_EQ(useHI, false);
	EXPECT_BOOL_EQ(useHII,false);
	EXPECT_BOOL_EQ(useH2, false);
	
	//get type density is independent from type activation, get density is not independent 
	//check if get density returns 0. (should give a error massage in log file)
	Hges = n.getDensity(p);
	EXPECT_DOUBLE_EQ(Hges, 0);
}


	

} //namespace crpropa
