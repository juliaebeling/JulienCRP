#include "crpropa/Massdistribution/Massdistribution.h"
#include "crpropa/Massdistribution/Cordes.h"
#include "crpropa/Massdistribution/Ferriere07.h"
#include "crpropa/Massdistribution/Nakanishi.h"
#include "crpropa/Massdistribution/Pohl2008.h"
#include "crpropa/Massdistribution/ConstantDensity.h"
#include "crpropa/Units.h"

#include "gtest/gtest.h"
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
	//check if all types are activited
	EXPECT_TRUE(useHI); 
	EXPECT_TRUE(useHII);
	EXPECT_TRUE(useH2);
	
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
	EXPECT_FALSE(useHI);
	EXPECT_FALSE(useHII);
	EXPECT_FALSE(useH2);
	
	//get type density is independent from type activation, get density is not independent 
	//check if get density returns 0. (should give a error massage in log file)
	Hges = n.getDensity(p);
	EXPECT_DOUBLE_EQ(Hges, 0);
}

TEST(testMassdistribution, SimpleTest) {
	
	Massdistribution MD;
	MD.add(new constantDensity(1,1,1));	//only HI types
	MD.add(new constantDensity(0,2,0));	//only HII type
	MD.add(new constantDensity(0,0,3));	//only H2 type

	Vector3d p(50*pc, 20*pc, -100*pc);	//random position for testing density 
	
	//check get density output
	double HI = MD.getHIDensity(p);
	double HII= MD.getHIIDensity(p);
	double H2 = MD.getH2Density(p);
	double n = MD.getDensity(p);
	EXPECT_DOUBLE_EQ(HI,1);
	EXPECT_DOUBLE_EQ(HII,2);
	EXPECT_DOUBLE_EQ(H2,3);
	EXPECT_DOUBLE_EQ(n,6);	//total density 2+2+3
	
	//check get function for type 
	EXPECT_TRUE(MD.getisforHI());
	EXPECT_TRUE(MD.getisforHII());
	EXPECT_TRUE(MD.getisforH2());

} 

TEST(testMassdistributionSuperposition, SimpleTest) {

	MassdistributionSuperposition MS;
	MS.addDensity(new constantDensity(1,1,2));	//sum 4
	MS.addDensity(new constantDensity(2,3,1));	//sum 6
	
	Vector3d p(50*pc,10*pc,-30*pc);	//random position for testing density
	EXPECT_DOUBLE_EQ(MS.getHIDensity(p),3);
	EXPECT_DOUBLE_EQ(MS.getHIIDensity(p),4);
	EXPECT_DOUBLE_EQ(MS.getH2Density(p),3);
	EXPECT_DOUBLE_EQ(MS.getDensity(p),10);	//sum of sums
	
}

TEST(testCordes, SimpleTest) {

	Cordes n;
	
	//check type Information
	EXPECT_FALSE(n.getisforHI());
	EXPECT_TRUE(n.getisforHII());
	EXPECT_FALSE(n.getisforH2());

	Vector3d p(3.1*kpc,2.9*kpc,-30*pc);	//position for testing density	
	
	double HII = n.getHIIDensity(p);
	double Hges = n.getDensity(p);
	
	EXPECT_NEAR(HII, 184500.,1);	// output in m^-3 ; uncertainty of 1e-6 cm^-1 
	EXPECT_NEAR(Hges, 184500.,1);	
	p.z=30*pc;			// invariant density for +/- z
	EXPECT_DOUBLE_EQ(HII,n.getDensity(p)); 
}

TEST(testNakanishi, SimpleTest) {

	Nakanishi n;
	
	//check type Information
	EXPECT_TRUE(n.getisforHI());
	EXPECT_FALSE(n.getisforHII());
	EXPECT_TRUE(n.getisforH2());
	
	Vector3d p(4*kpc,-2.5*kpc,-0.85*kpc);	//position for testing density
	
	//testing HI component
	EXPECT_NEAR(n.getHIPlanedensity(p),162597,1);	//uncertaincy of 1e-6 cm^-3 
	EXPECT_NEAR(n.getHIScaleheight(p),0.3109*kpc,0.1*pc); 
	EXPECT_NEAR(n.getHIDensity(p),914,1); 	//uncertainc 1e-6 cm^-3 
	
	//testing HII compontent
	EXPECT_NEAR(n.getH2Planedensity(p),741999,1); //uncertaincy of 1e-6 cm^-3
	EXPECT_NEAR(n.getH2Scaleheight(p),88.2*pc,0.1*pc);
	EXPECT_NEAR(n.getH2Density(p),0,1);
	
	//testing total Density
	EXPECT_NEAR(n.getDensity(p),914,2); //double uncertaincy for both type á 1cm^-3
	
	
	
	p = Vector3d(50*pc,100*pc,10*pc);	// second position for testing density
	
	//testing HI component
	EXPECT_NEAR(n.getHIPlanedensity(p),543249,1);
	EXPECT_NEAR(n.getHIScaleheight(p),125.6*pc,0.1*pc);
	EXPECT_NEAR(n.getHIDensity(p),540867,1);
	
	//testing H2 component
	EXPECT_NEAR(n.getH2Planedensity(p),10556748,1);
	EXPECT_NEAR(n.getH2Scaleheight(p),57.2*pc,0.1*pc);
	EXPECT_NEAR(n.getH2Density(p),10335137,1);
	
	//test set type function
	n.setisforHI(false);
	EXPECT_FALSE(n.getisforHI());
	EXPECT_TRUE(n.getisforH2());
	
	n.setisforH2(false);
	EXPECT_FALSE(n.getisforHI());
	EXPECT_FALSE(n.getisforH2());
	
	//check if density output is zero if all density-types are deaktivated (should give warning in log-file)
	EXPECT_DOUBLE_EQ(n.getDensity(p),0);	
}

TEST(testFerriere, SimpleTest) {
	
	Ferriere n;
	
	//check type information
	
	EXPECT_TRUE(n.getisforHI());
	EXPECT_TRUE(n.getisforHII());
	EXPECT_TRUE(n.getisforH2());
	
	//testing density in inner Ring (R <= 3*kpc)
	Vector3d p(-60*pc,60*pc,-20*pc);	//testing position in region of CMZ
	
	//test CMZ Trafo
	Vector3d Trafo;
	Trafo = n.CMZTrafo(p);
	EXPECT_NEAR(Trafo.x,5.9767*pc,1e-4*pc);
	EXPECT_NEAR(Trafo.y,12.8171*pc,1e-4*pc);
	EXPECT_DOUBLE_EQ(Trafo.z,p.z);	//no transformation in z component
	
	//test DISK Trafo
	Trafo = n.DISKTrafo(p);
	EXPECT_NEAR(Trafo.x,11.0660*pc,1e-4*pc);
	EXPECT_NEAR(Trafo.y,82.5860*pc,1e-4*pc);
	EXPECT_NEAR(Trafo.z,-25.6338*pc,1e-4*pc);
	
	//testing density 
	EXPECT_NEAR(n.getHIDensity(p),6237723,1); 	//uncertaincy 1e-6 cm^-3
	EXPECT_NEAR(n.getH2Density(p),35484825,1); 
	EXPECT_NEAR(n.getHIIDensity(p),5570335,1);	
	EXPECT_NEAR(n.getDensity(p),47292883,1);
	
	Vector3d p2(500*pc,900*pc,35*pc);	//testing position in region of the DISK
	EXPECT_NEAR(n.getHIIDensity(p2),48190,1);
	EXPECT_NEAR(n.getHIDensity(p2),5,1);
	EXPECT_NEAR(n.getH2Density(p2),0,1);
	EXPECT_NEAR(n.getDensity(p2),48195,1);
	
	//testing the outer region R>3kpc
	
	Vector3d p3(5*kpc,4*kpc,-29*pc);	//testing position with 3kpc < R < R_sun
	EXPECT_NEAR(n.getHIDensity(p3),540607,1);
	EXPECT_NEAR(n.getHIIDensity(p3),66495 ,1);
	EXPECT_NEAR(n.getH2Density(p3),2492685,1);
	EXPECT_NEAR(n.getDensity(p3),3099787,1);
	
	Vector3d p4(10*kpc,2*kpc,50*pc);	//testing position with R > R_sun
	EXPECT_NEAR(n.getHIDensity(p4),431294,1);
	EXPECT_NEAR(n.getHIIDensity(p4),22109,1);
	EXPECT_NEAR(n.getH2Density(p4),54099,1);
	EXPECT_NEAR(n.getDensity(p4),507502,1);
	
	//test get/set type funktion
	
	n.setisforHI(false);
	EXPECT_FALSE(n.getisforHI());
	EXPECT_TRUE(n.getisforHII());
	EXPECT_TRUE(n.getisforH2());
	
	n.setisforHII(false);
	EXPECT_FALSE(n.getisforHI());
	EXPECT_FALSE(n.getisforHII());
	EXPECT_TRUE(n.getisforH2());
	
	n.setisforH2(false);
	EXPECT_FALSE(n.getisforHI());
	EXPECT_FALSE(n.getisforHII());
	EXPECT_FALSE(n.getisforH2());
	
	//check if density is set to zero if all types are deaktivated (schould give warning in log-file)
	EXPECT_DOUBLE_EQ(n.getDensity(p),0);
}

TEST(testPohl,SimpleTest) {
	
	Pohl08 n;
	
	n.loadGridHI();
	n.loadGridH2();

	//check type information
	EXPECT_TRUE(n.getisforHI());
	EXPECT_FALSE(n.getisforHII());
	EXPECT_TRUE(n.getisforH2());

	//test boundary of grid
	Vector3d px(50*kpc,10*kpc,5*pc);
	EXPECT_DOUBLE_EQ(n.getHIDensity(px),0);
	EXPECT_DOUBLE_EQ(n.getH2Density(px),0);
	
	Vector3d py(5*kpc,30*kpc,0.5*kpc);
	EXPECT_DOUBLE_EQ(n.getHIDensity(py),0);
	EXPECT_DOUBLE_EQ(n.getH2Density(py),0);
	
	//test set type funktion
	n.setisforHI(false);
	EXPECT_FALSE(n.getisforHI());
	EXPECT_FALSE(n.getisforHII());
	EXPECT_TRUE(n.getisforH2());
	
	n.setisforH2(false);
	EXPECT_FALSE(n.getisforHI());
	EXPECT_FALSE(n.getisforHII());
	EXPECT_FALSE(n.getisforH2());
	
	//check if density is set to zero when all densties are deaktivated (sould give Warning in log-File)
	Vector3d p(0.);
	EXPECT_DOUBLE_EQ(n.getDensity(p),0);
	
}

} //namespace crpropa
