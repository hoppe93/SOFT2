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
const unsigned int 
    IGAMMA = 1,
    IZEFF = 0;

const unsigned int Test_BremsstrahlungEmission::NTESTVALUES = 4;

const slibreal_t Test_BremsstrahlungEmission::TEST_ZEFF_GAMMA[4][2] = {
{2,	    5},
{7,     13},
{10,    38},
{15,	67}
};

const slibreal_t Test_BremsstrahlungEmission::TEST_VAL[NTESTVALUES][NTESTVALUES] = {
{2.88228628200495e-30,	9.39765963487424e-30,	2.10302400717392e-29,	2.84214668681879e-29},
{3.53080069545606e-29,	1.15121330527209e-28,	2.57620440878806e-28,	3.48162969135301e-28},
{7.20571570501237e-29,	2.34941490871856e-28,	5.25756001793481e-28,	7.10536671704696e-28},
{1.62128603362778e-28,	5.28618354461676e-28,	1.18295100403533e-27,	1.59870751133557e-27}
};


/**
 * Check that the formula for total emission
 * of synchrotron radiation has been implemented
 * correctly. Tests vs values obtained from same formulas implemented in matlab
 */
/*bool Test_BremsstrahlungEmission::CheckTotalEmission(const slibreal_t tol) {
    unsigned int i;
    Detector *det = GetDetector(0);
    slibreal_t pwr, corr, Delta;

    for (i = 0; i < NTESTVALUES; i++) {
	    slibreal_t Input_Z[6] = {TEST_Z[i][0], TEST_Z[i][1], TEST_Z[i][2], TEST_Z[i][3], TEST_Z[i][4], TEST_Z[i][5]};
        slibreal_t Input_DENS[6] = {TEST_DENS[i][0], TEST_DENS[i][1], TEST_DENS[i][2], TEST_DENS[i][3], TEST_DENS[i][4], TEST_DENS[i][5]};
        slibreal_t Input_Z0[6] = {TEST_Z0[i][0], TEST_Z0[i][1], TEST_Z0[i][2], TEST_Z0[i][3], TEST_Z0[i][4], TEST_Z0[i][5]};
        ConeBremsstrahlungEmission cbse(det, nullptr, TEST_NR_POW[i].nspecies, Input_Z, Input_Z0, Input_DENS);

        cbse.CalculateTotalEmission();
        pwr = cbse.GetTotalEmission();
        corr = TEST_NR_POW[i].val;
        
        Delta = fabs((pwr-corr)/corr);

        if (Delta >= tol) {
            this->PrintError("Total bremsstrahlung  emission was not calculated correctly. Delta = %e", Delta);
            return false;
        }
    }
    return true;
}*/

/* 
 * Test so that the spectrum integrates to the correct value
 */
bool Test_BremsstrahlungEmission::CheckSpectrumEmission(const slibreal_t tol) {
    
    unsigned int i, j;
    Detector *det = GetDetector(50, 1, 50); //Gives detector and defines which photon momenta to look at
    MagneticFieldAnalytical2D *dummy_mf = GetMagneticField();
    slibreal_t pwr, corr, Delta;

    /*slibreal_t *Z = new slibreal_t[nspecies],
    *Z0 = new slibreal_t[nspecies], 
    *DENS = new slibreal_t[nspecies];
    for (i = 0; i < nspecies; i++){
        Z[i] = TEST_Z[NTESTVALUES-1][i];
        Z0[i] = TEST_Z0[NTESTVALUES-1][i];
        DENS[i] = TEST_DENS[NTESTVALUES-1][i];
    }*/
    for (j = 0; j < NTESTVALUES; j++) {
        ConeBremsstrahlungEmission cbe(det, nullptr, TEST_ZEFF_GAMMA[j][IZEFF]);
        for (i = 0; i < NTESTVALUES; i++) {
            
            RadiationParticle *rp = GetRadiationParticle(i, det, dummy_mf);
            
            cbe.CalculateSpectrum(rp);
            cbe.IntegrateSpectrum();
            pwr = cbe.GetTotalEmission();
            corr = TEST_VAL[j][i];

            
            Delta = fabs((pwr-corr)/corr);

            if (Delta >= tol) {
                this->PrintError("Total bremsstrahlung emission was not calculated correctly. Delta = %e", Delta);
                return false;
            }
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
        visang   = 1.0,
        dir[3] = {0.0,1.0,0.0},
        pos[3] = {0.0,-1.069,0.0};
        
    Vector<3>
        direction(dir),
        position(pos);

    return new Detector(aperture, visang, direction, position, nwavelengths, l0, l1);
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
        gamma = TEST_ZEFF_GAMMA[i][IGAMMA],
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
        this->PrintOK("Total emission is implemented correctly.");
    else
        success = false;*/

    if (CheckSpectrumEmission(SPECT_TOL))
        this->PrintOK("Spectrum integrates correctly.");
    else
        success = false;
    
    return success;
}

