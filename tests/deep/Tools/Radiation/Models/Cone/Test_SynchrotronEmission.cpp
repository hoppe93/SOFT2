/**
 * Test of the synchrotron emission implemented in the cone model.
 * The following tests are included:
 *
 *   [ ] - Emission formulas integrate to Larmor's formula.
 *
 */

#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/Cone/ConeSynchrotronEmission.h"
#include "Test_SynchrotronEmission.h"

using namespace std;
using namespace __Radiation;

const unsigned int
    IGAMMA = 0,
    ICOORDINATES = 1;

const unsigned int Test_SynchrotronEmission::NTESTPARTICLES = 15;
const slibreal_t Test_SynchrotronEmission::TESTPARTICLES[NTESTPARTICLES][4] = {
    // gamma,  x,  	y,	z
    { 2.0,     0.10,    0.50,	0.00},
    { 2.0,     0.40,    0.20,	0.10},
    { 2.0,     0.90,    0.00,	0.00},

    { 5.0,     0.15,    0.60, 	0.20},
    { 5.0,     0.45,    0.50,	-0.20},
    { 5.0,     0.70,    0.30,	-0.20},

    {15.0,     -0.20,    0.50,	0.00},
    {15.0,     -0.50,    0.00,	0.10},
    {15.0,     0.80,    -0.10,	-0.05},

    {55.0,     0.40,    0.40,	0.00},
    {55.0,     0.55,    0.10,	0.10},
    {55.0,     -0.75,    0.00,	0.20},

    {100.0,    -0.00,   -0.68,	0.00},
    {100.0,    0.00,    0.50,	-0.10},
    {100.0,    0.00,    0.85,	0.00}
};

/**
 * Larmor's formula.
 */
slibreal_t Test_SynchrotronEmission::Larmor(
    RadiationParticle *rp
) {
    const slibreal_t
        c = LIGHTSPEED,
        e = ELEMENTARY_CHARGE,
        eps0 = EPS0,
        m = ELECTRON_MASS;

    slibreal_t
        ppar = rp->GetPpar(),
        pperp = rp->GetPperp(),
        p2 = rp->GetP2(),
        p = sqrt(p2),
        gamma = rp->GetGamma(),
        gammapar = gamma / sqrt(1.0 + pperp*pperp),
        beta = p / gamma,
        cosThetap = ppar / p,
        betaperp = pperp / gamma,
        B = rp->GetB();
    
    slibreal_t C = e*e*e*e / (6.0*M_PI*eps0*m*m*c);
    slibreal_t Cp = B*B * gamma*gamma * gammapar*gammapar * betaperp*betaperp * (1.0 - beta*cosThetap*cosThetap);

    return C * Cp;
    /*slibreal_t C = e*e * B*B * gamma*gamma*gamma*gamma * gammapar*gammapar * betaperp*betaperp /
        (2.0*sqrt(3.0)*M_PI*m*m*c*c);

    return C;*/
}

/**
 * Check that the formula for total emission
 * of synchrotron radiation has been implemented
 * correctly.
 */
bool Test_SynchrotronEmission::CheckTotalEmission(const slibreal_t tol) {
    unsigned int i;
    Detector *det = GetDetector(0);
    ConeSynchrotronEmission cse(det, nullptr);
    slibreal_t pwr, corr, Delta;
	MagneticFieldAnalytical2D *mf = GetMagneticField();

    for (i = 0; i < NTESTPARTICLES; i++) {
        RadiationParticle *rp = GetRadiationParticle(i, det, mf);
        cse.CalculateTotalEmission(rp);
        pwr = cse.GetTotalEmission();
        corr = Larmor(rp);
        
        Delta = fabs((pwr-corr)/corr);
        if (Delta >= tol) {
            this->PrintError("Total synchrotron emission was not calculated correctly. Delta = %e", Delta);
            return false;
        }
    }

    return true;
}

/**
 * Check that the integrated spectrum has been
 * implemented correctly.
 */
bool Test_SynchrotronEmission::CheckSpectrumEmission(const slibreal_t tol) {
    unsigned int i;
    slibreal_t pwr, corr, Delta, lambdac;
    Detector *dummyDet = GetDetector(0);

	MagneticFieldAnalytical2D *dummy_mf = GetMagneticField();

    for (i = 0; i < NTESTPARTICLES; i++) {
        RadiationParticle *rp = GetRadiationParticle(i, dummyDet, dummy_mf); 

        lambdac = GetLambdaC(rp);

        Detector *det = GetDetector(40001, lambdac/10.0, lambdac*1000.0);
        ConeSynchrotronEmission cse(det, dummy_mf);

        cse.CalculateSpectrum(rp);
        cse.IntegrateSpectrum();
        
        pwr = cse.GetTotalEmission();
        corr = Larmor(rp);

        Delta = fabs((pwr-corr)/corr);
        if (Delta >= tol) {
            this->PrintError("Total synchrotron emission was not calculated correctly for test particle #%u. Delta = %e", i, Delta);
            return false;
        }
    }

    return true;
}

/**
 * Returns a sample detector.
 */
Detector *Test_SynchrotronEmission::GetDetector(unsigned int nwavelengths, slibreal_t l0, slibreal_t l1) {
    const slibreal_t
        aperture = 0.006,
        visang   = 1.0,
        dir[3] = {0.0,1.0,0.0},
        pos[3] = {0.0,-1.069,0.0};
        
    Vector<3>
        direction(dir),
        position(pos);

    return new Detector(aperture, visang, direction, position, nwavelengths, l0, l1);
}

/**
 * Returns the critical wavelength for the synchrotron spectrum
 * emitted by the given particle.
 */
slibreal_t Test_SynchrotronEmission::GetLambdaC(RadiationParticle *rp) {
    slibreal_t
        gamma = rp->GetGamma(),
        pperp = rp->GetPperp(),
        gamma2 = gamma*gamma,
        ppar2 = rp->GetPpar() * rp->GetPpar(),
        pperp2 = pperp*pperp,
        betapar2 = ppar2 / gamma2,
        betaperp2 = pperp2 / gamma2,
        beta = sqrt(betapar2+betaperp2),
        l = 4.0*M_PI*ELECTRON_MASS*beta*LIGHTSPEED / (3.0*gamma2*ELEMENTARY_CHARGE);
    	
        Vector<3>& rhat = rp->GetRCP();
	    rhat.Normalize();
	    Vector<3> pos = rp->GetPosition();
        Vector<3> Bvec = rp->GetBvec();
	    Vector<3> rh_cr_rh_cr_Bvec = Vector<3>::Cross(rhat, Vector<3>::Cross(rhat, Bvec));

        slibreal_t lc = l*1/rh_cr_rh_cr_Bvec.Norm();

    return lc;
}

/**
 * Returns the pre-defined RadiationParticle object
 * with index 'i'.
 *
 * i: Index in table 'predef_particles' with pre-defined
 *    particle parameters to test.
 */
RadiationParticle *Test_SynchrotronEmission::GetRadiationParticle(unsigned int i, Detector *det, MagneticFieldAnalytical2D *mf) {
    if (i >= NTESTPARTICLES)
        throw SOFTException("Trying to access non-existant test-particle.");

    Vector<3> x (TESTPARTICLES[i]+ICOORDINATES), 
		rhat = x-det->GetPosition(),
		Bvec = mf->Eval(x[0], x[1], x[2]),
		bHat = Bvec;

	bHat.Normalize();
	rhat.Normalize();
		

    slibreal_t
        Jdtdrho = 1.0,
        gamma = TESTPARTICLES[i][IGAMMA],
        p2 = gamma*gamma - 1.0,
        p  = sqrt(p2),
        cosThetap = rhat.Dot(bHat),
        sinThetap = sqrt(1.0 - cosThetap*cosThetap),
        ppar = p * cosThetap,
        pperp = p * sinThetap,
        B = Bvec.Norm(),
        m = ELECTRON_MASS, q = -ELEMENTARY_CHARGE;
		

		
    Vector<3> P = p*rhat;

    return new RadiationParticle(
        x, P,
        Jdtdrho, ppar, pperp,
        gamma, p2, det->GetPosition(),
        B, Bvec, bHat, m, q, 0, 0, 0
    );
}

/**
 * Run all the tests in the module.
 */
bool Test_SynchrotronEmission::Run(bool) {
    bool success = true;
    slibreal_t TOTEM_TOL=SQRT_REAL_EPSILON,
               SPECT_TOL=1e-3;

    if (CheckTotalEmission(TOTEM_TOL))
        this->PrintOK("Total emission is implemented correctly.");
    else
        success = false;

    if (CheckSpectrumEmission(SPECT_TOL))
        this->PrintOK("Spectrum integrates to Larmor's formula.");
    else
        success = false;
    
    return success;
}

