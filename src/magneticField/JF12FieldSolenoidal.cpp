#include "crpropa/magneticField/JF12FieldSolenoidal.h"
#include "crpropa/Units.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"

namespace crpropa {

JF12FieldSolenoidal::JF12FieldSolenoidal(double delta, double zs) {
	zS = zs; // set scale heigth for the parabolic X field lines
	r1 = 5 * kpc; // inner boundary of the disk field
	r2 = 20 * kpc; // outer boudary of the disk field
	r1s = r1 + delta; // the magnetic flux of the spirals is redirected for r in [r1,r1s]
	r2s = r2 - delta; // same here at outer boundary between r2s and r2
	phi0 = 0.; // somewhat arbitrary choice, has to be chosen in [-pi,pi]

	for (int i = 1;i < 9; i++){
		// fill the array with angles in [-pi,pi] where the 8 spiral arms intersect the r1 - ring
		// indexing starts at 1 to match the indexing in the papers on the JF12 field!
		phi0Arms[i] = M_PI - cotPitch * log(rArms[i-1] / r1);
	}

	// cyclic closure of the array, with next values periodically continued
	// outside [-pi,pi] to simplify looping and searching for correct spiral arms
	phi0Arms[0] = phi0Arms[8] + 2 * M_PI;
	phi0Arms[9] = phi0Arms[1] - 2 * M_PI;
	phi0Arms[10] = phi0Arms[2] - 2 *M_PI;

	// determine the position of phi0 in the array, i.e. find the correct spiral arm.
	int idx0 = 1; // corresponding index in phi0Arms such that phi0Arms[idx0] < phi0 < phi0Arms[idx0-1]
	while (phi0 < phi0Arms[idx0]){
		idx0 += 1; // search clockwise, starting with the check if phi0Arms[1] < phi0 < phi0Arms[0]
	}

	// fill the bDisk array with spiral field strengths at r = r1.
	// note the indexing starting with 1 here to match the indexing in the JF12 papers!
	// for a position (r1,phi), phi in [-pi,pi], the correct field strength is given by
	// bDisk[i] if phi0Arms[i] < phi0 < phi0Arms[i-1].
	bDisk[1] = 0.1 * muG;
	bDisk[2] = 3.0 * muG;
	bDisk[3] = -0.9 * muG;
	bDisk[4] = -0.8 * muG;
	bDisk[5] = -2.0 * muG;
	bDisk[6] = -4.2 * muG;
	bDisk[7] = 0.0 * muG;

	// re-compute b_8 for actual (net flux = 0)-correction of the spiral field with minimal round-off errors
	double flux1to7 = 0.;
	for (int i = 1; i < 8; i++){
		flux1to7 += (phi0Arms[i-1] - phi0Arms[i]) * bDisk[i];
	}
	bDisk[8] = -flux1to7 / (phi0Arms[7] - phi0Arms[8]);

	bDisk[0] = bDisk[8]; // again close the array periodically
	bDisk[9] = bDisk[1];
	bDisk[10] = bDisk[2];

	// set coefficients for the evaluation of the phi-integral over the piecewise constant field strengths at r=r1
	// such that it may be evaluated as H(phi) = phiCoeff[j] + bDisk[j] * phi later on
	// start integration at phi0Arms[0] first, shift to lower integration boundary phi0 later
	phiCoeff[0] = 0;
	for (int i = 1; i < 10; i++){
		phiCoeff[i] = phiCoeff[i-1] + (bDisk[i-1] - bDisk[i]) * phi0Arms[i-1];
	}

	// correct for H(phi0) = 0
	corr = phiCoeff[idx0] + bDisk[idx0] * phi0;
	for (int i = 1; i < 10; i++){
		phiCoeff[i] = phiCoeff[i] - corr;
	}
}

void JF12FieldSolenoidal::setDiskTransitionWidth(double delta) {
	r1s = r1 + delta;
	r2s = r2 - delta;
}

void JF12FieldSolenoidal::setXScaleHeight(double zs) {
	zS = zs;
}

double JF12FieldSolenoidal::getDiskTransitionWidth() const {
	return (r1s - r1);
}

double JF12FieldSolenoidal::getXScaleHeight() const {
	return zS;
}

void JF12FieldSolenoidal::deactivateOuterTransition() {
	r2s = r2;
}

void JF12FieldSolenoidal::setUseStriatedField(bool use) {
	if ((use) and (striatedGrid)) {
		KISS_LOG_WARNING << "JF12FieldSolenoidal: No striated field set: ignored.";
		return;
	}
	useStriatedField = use;
}

void JF12FieldSolenoidal::setUseTurbulentField(bool use) {
	if ((use) and (turbulentGrid)) {
		KISS_LOG_WARNING << "JF12FieldSolenoidal: No turbulent field set: ignored.";
		return;
	}
	useTurbulentField = use;
}

Vector3d JF12FieldSolenoidal::getDiskField(const double& r, const double& z, const double& phi, const double& sinPhi, const double& cosPhi) const {
	Vector3d b(0.);

	if (useDiskField){
		double lfDisk = logisticFunction(z, hDisk, wDisk); // for vertical scaling as in initial JF12

		double hint = getHPhiIntegral(r, phi); // phi integral to restore solenoidality in transition region, only enters if r is in [r1,r1s] or [r2s,r2]
		double mag1 = getSpiralFieldStrengthConstant(r, phi); // returns bDisk[j] for the current spiral arm

		if ((r1 < r) && (r < r2)) {
			double pdelta = getDiskTransitionPolynomial(r);
			double qdelta = getDiskTransitionPolynomialDerivative(r);
			double br = pdelta * mag1 * sinPitch;
			double bphi = pdelta * mag1 * cosPitch - qdelta * hint * sinPitch;

			b.x += br * cosPhi - bphi * sinPhi;
			b.y += br * sinPhi + bphi * cosPhi;

			b *= (1 - lfDisk);
		}
	}
	return b;
}

Vector3d JF12FieldSolenoidal::getXField(const double& r, const double& z, const double& sinPhi, const double& cosPhi) const {
	Vector3d b(0.);

	if (useXField){
		double bMagX;
		double sinThetaX, cosThetaX;
		double rp; // radius where current intial field line passes z = 0
		double rc = rXc + fabs(z) / tanThetaX0;
		double r0c = rXc + zS / tanThetaX0; // radius where field line through rXc passes z = zS
		double f, r0, br0, bz0;
		bool inner = true; // distinguish between inner and outer region

		// return intial field if z>=zS
		if (fabs(z) > zS){
			if ((r == 0.)){
				b.z = bX / ((1. + fabs(z) * cotThetaX0 / rXc) * (1. + fabs(z) * cotThetaX0 / rXc));
				return b;
			}

			if (r < rc) {
			// inner varying elevation region
				rp = r * rXc / rc;
				bMagX = bX * exp(-1 * rp / rX) * (rXc / rc) * (rXc / rc);

				double thetaX = atan(fabs(z) / (r - rp));

				if (z == 0)
					thetaX = M_PI / 2.;

				sinThetaX = sin(thetaX);
				cosThetaX = cos(thetaX);
			}
			else {
			// outer constant elevation region
				rp = r - fabs(z) / tanThetaX0;
				bMagX = bX * exp(-rp / rX) * (rp / r);

				sinThetaX = sinThetaX0;
				cosThetaX = cosThetaX0;
			}
			double zsign = z < 0 ? -1 : 1;
			b.x += zsign * bMagX * cosThetaX * cosPhi;
			b.y += zsign * bMagX * cosThetaX * sinPhi;
			b.z += bMagX * sinThetaX;
		}
		// parabolic field lines for z<zS
		else {
				// determine r at which parabolic field line through (r,z) passes z = zS
				r0 = r * 1. / (1.- 1./ (2. * (zS + rXc * tanThetaX0)) * (zS - z * z / zS));

				// determine correct region (inner/outer)
				if (r0 >= r0c){
					r0 = r + 1. / (2. * tanThetaX0) * (zS - z * z / zS);
					inner = false;
				}

				// field strength at that position
				 if (r0 < r0c){
					 rp = r0 * rXc / r0c;
					 double thetaX = atan(zS / (r0 - rp));

					 // field strength at (r0,zS) for inner region
					 br0 = bX * exp(- rp / rX) * (rXc/ r0c) * (rXc/ r0c) * cos(thetaX);
					 bz0 = bX * exp(- rp / rX) * (rXc/ r0c) * (rXc/ r0c) * sin(thetaX);
				 }
				 else {
					 // field strength at (r0,zS) for outer region
					 rp = r0 - zS / tanThetaX0;
					 br0 =  bX * exp(- rp / rX) * (rp/r0) * cosThetaX0;
					 bz0 =  bX * exp(- rp / rX) * (rp/r0) * sinThetaX0;
				 }

				 // compute factor F for solenoidality
				 if (inner){
					 f = 1. / ((1. - 1./( 2. + 2. * (rXc * tanThetaX0/ zS)) * (1. - (z / zS) * (z / zS))) * (1. - 1./( 2. + 2. * (rXc * tanThetaX0/ zS)) * (1. - (z / zS) * (z / zS))));
				 }
				 else {
					 f = 1. + 1/ (2 * r * tanThetaX0/ zS) * (1. - (z / zS) * (z / zS));
				 }

				 double br = z / zS * f * br0;
				 double bz = bz0 * f;

				 b.x += br * cosPhi;
				 b.y += br * sinPhi;
				 b.z += bz;
		}
	}
	return b;
}

double JF12FieldSolenoidal::getDiskTransitionPolynomial(const double& r) const {
	// 0 disk field outside
	if ((r < r1) || (r > r2)) {
		return 0.;
	}
	// unchanged field
	if ((r > r1s) && (r < r2s)) {
		return r1/r;
	}
	// transitions region parameters
	double r_a = r1;
	double r_b = r1s;

	if (r >= r2s) {
		r_a = r2;
		r_b = r2s;
	}
	// differentiable transition at r_s, continous at r_a
	double fakt = (r_a / r_b - 2.) / ((r_a - r_b) *  (r_a - r_b));
	return (r1/r_b) * (2. - r / r_b + fakt * (r-r_b) * (r-r_b));
}

double JF12FieldSolenoidal::getDiskTransitionPolynomialDerivative(const double& r) const {
	// 0 disk field outside
	if ((r < r1) || (r > r2)) {
		return 0.;
	}
	// unchanged field
	if ((r > r1s) && (r < r2s)) {
		return 0.;
	}
	// transitions region parameters
	double r_a = r1;
	double r_b = r1s;

	if (r >= r2s) {
		r_a = r2;
		r_b = r2s;
	}
	// differentiable transition polynomial at r_s, continous at r_a
	double fakt = (r_a / r_b - 2.) / ((r_a - r_b) * (r_a - r_b));
	return (r1/r_b) * (2. - 2. * r/r_b + fakt * (3. * r * r - 4. * r * r_b + r_b * r_b));
}

double JF12FieldSolenoidal::getHPhiIntegral(const double& r, const double& phi) const {
	// Evaluates the H(phi1) integral for solenoidality for the position (r,phi) which is mapped back to (r1=5kpc,phi1)
	// along the spiral field line.
	double H_ret = 0.;
	int idx = 1;

	if ((r1 < r) && (r < r2)){
		// find index of the correct spiral arm for (r1,phi1) just like in getSpiralFieldStrengthConstant
		double phi1 = phi - log(r/r1) * cotPitch;
		phi1 = atan2(sin(phi1), cos(phi1));
		while (phi1 < phi0Arms[idx]){
			idx += 1;
		}
		H_ret = phi1 * bDisk[idx] + phiCoeff[idx];
	}
	return H_ret;
}

double JF12FieldSolenoidal::getSpiralFieldStrengthConstant(const double& r, const double& phi) const {
	// For a given position (r, phi) in polar coordinates, this method returns the field strength
	// of the spiral field at r1 = 5 kpc for the magnetic spiral arm where (r, phi) is located.
	// The method first computes the angle phi1 at which the spiral field line passing through (r, phi) intersects
	// the circle with radius r1 = 5 kpc. Afterwards, the correct spiral arm is found by searching the index idx
	// such that phi0Arms[idx] < phi1 < phi0Arms[idx-1]. The correct field strength of the respective spiral arm
	// where (r, phi) is located is then given as bDisk[idx].
	double b_ret = 0.;
	int idx = 1;
	if ((r1 < r) && (r < r2)){
		double phi1 = phi - log(r/r1) * cotPitch; // map the position (r, phi) to (5 kpc, phi1) along the logarithmic spiral field line
		phi1 = atan2(sin(phi1), cos(phi1)); // map this angle to [-pi,+pi]
		while (phi1 < phi0Arms[idx]){
			idx += 1; // run clockwise through the spiral arms; the cyclic closure of phi0Arms[9] = phi0Arms[1] - 2 pi is needed if -pi <= phi1 <= phi0Arms[8].
		}
		b_ret = bDisk[idx];
	}
	return b_ret;
}
} // namespace crpropa
