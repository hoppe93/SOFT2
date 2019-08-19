/**
 * Test of the bremsstrahlung screned emission implemented in the cone model.
 * The following tests are included:
 *
 *   [ ] - Emission formulas (Only integrated spectrum).  
 *
 */

#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/Cone/ConeBremsstrahlungEmission.h"
#include "Test_BremsstrahlungEmission.h"

using namespace std;
using namespace __Radiation;

//Constants used to retrive data from tables
const unsigned int IZ = 0,
    IDENS = 1,
    IGAMMA = 2,
    IINTVAL = 3;

const unsigned int Test_BremsstrahlungEmission::NTESTVALUES = 4;

const slibreal_t Test_BremsstrahlungEmission::TEST_VALUES[NTESTVALUES][NTESTVALUES]{
{2,	    1.64900000000000e+19,	5,	    2.41757274275679e-08},
{7,	    1.14670000000000e+20,	13,	    7.88246674899129e-08},
{10,	7.39030000000000e+19,	38,	    1.76395160635127e-07},
{15,	9.10030000000000e+19,	67,	    2.38390488962464e-07}
};

/**
 * Check that the formula for total emission
 * of synchrotron radiation has been implemented
 * correctly. Tests vs values obtained from same formulas implemented in matlab
 */
/*bool Test_BremsstrahlungEmission::CheckTotalEmission(const slibreal_t tol) { //Needs to be implemented
   
}*/

/* 
 * Test so that the spectrum integrates to the correct value
 */
bool Test_BremsstrahlungEmission::CheckSpectrumEmission(const slibreal_t tol) {
    
    Detector *det = GetDetector(50, 1, 50); //Gives detector and defines which photon momenta to look at
    MagneticFieldAnalytical2D *dummy_mf = GetMagneticField();
    slibreal_t pwr, corr, Delta;

    unsigned int nspecies = NTESTVALUES;
    slibreal_t *Z = new slibreal_t[nspecies], 
        *dens = new slibreal_t[nspecies];

    for (unsigned int i = 0; i < nspecies; i++){
        Z[i] = TEST_VALUES[i][IZ];
        dens[i] = TEST_VALUES[i][IDENS];
    }

    ConeBremsstrahlungEmission cbe(det, nullptr, nspecies, Z, dens);
    for (unsigned int i = 0; i < NTESTVALUES; i++) {
            
        RadiationParticle *rp = GetRadiationParticle(i, det, dummy_mf);
            
        cbe.CalculateSpectrum(rp);
        cbe.IntegrateSpectrum();
        pwr = cbe.GetTotalEmission();
        corr = TEST_VALUES[i][IINTVAL];
            
        Delta = fabs((pwr-corr)/corr);

        if (Delta >= tol) {
            this->PrintError("Total bremsstrahlung emission was not calculated correctly. Delta = %e", Delta);
            return false;
        }
    }
    return true;
}

/**
 * Returns a sample detector.
 */
Detector *Test_BremsstrahlungEmission::GetDetector(unsigned int nwavelengths, slibreal_t l0, slibreal_t l1) {
    const slibreal_t
        aperture = 0.006,
        roll     = 0,
        visang   = 1.0,
        dir[3]   = {0.0,1.0,0.0},
        pos[3]   = {0.0,-1.069,0.0};

        
    Vector<3>
        direction(dir),
        position(pos);

    return new Detector(aperture, roll, visang, direction, position, nwavelengths, l0, l1);
}

/**
 * Returns the pre-defined RadiationParticle object
 * with index 'i'.
 *
 * i: Index in table 'predef_particles' with pre-defined
 *    particle parameters to test, in this case gamma.
 */
RadiationParticle *Test_BremsstrahlungEmission::GetRadiationParticle(unsigned int i, Detector *det, MagneticFieldAnalytical2D *mf) {
    if (i >= NTESTVALUES)
        throw SOFTException("Trying to access non-existant test-value.");

    Vector<3> x, //Arbitrary position
		rhat = x-det->GetPosition(),
		Bvec = mf->Eval(x[0], x[1], x[2]),
		bHat = Bvec;

	bHat.Normalize();
	rhat.Normalize();
		

    slibreal_t
        Jdtdrho = 1.0,
        gamma = TEST_VALUES[i][IGAMMA],
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
bool Test_BremsstrahlungEmission::Run(bool) {
    bool success = true;
    //slibreal_t TOTEM_TOL=SQRT_REAL_EPSILON;
    slibreal_t SPECT_TOL=1e-3;

    /*if (CheckTotalEmission(TOTEM_TOL))
        this->PrintOK("Total emission is implemented correctly."); //For test of total emission
    else
        success = false;*/

    if (CheckSpectrumEmission(SPECT_TOL))
        this->PrintOK("Spectrum integrates correctly.");
    else
        success = false;
    
    return success;
}

