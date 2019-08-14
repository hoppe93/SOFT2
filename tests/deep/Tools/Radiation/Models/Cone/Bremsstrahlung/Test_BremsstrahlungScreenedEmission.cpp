/**
 * Test of the bremsstrahlung screned emission implemented in the cone model.
 * The following tests are included:
 *
 *   [ ] - Emission formulas.
 *
 */

#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/Cone/ConeBremsstrahlungScreenedEmission.h"
#include "Test_BremsstrahlungScreenedEmission.h"

using namespace std;
using namespace __Radiation;

//Constants used to retrive data from tables
const unsigned int 
    IGAMMA = 0,
    IINTSPEC = 1;

const unsigned int Test_BremsstrahlungScreenedEmission::NTESTVALUES = 10;

const slibreal_t Test_BremsstrahlungScreenedEmission::TEST_Z[NTESTVALUES][6] = {
{8}, {5}, {3, 2}, {5, 5}, {6, 2, 6, 5}, {2, 8, 5, 4}, {8, 6, 4, 7, 5}, {5, 6, 4, 3, 5}, {7, 4, 2, 4, 3, 4}, {7, 4, 4, 4, 2, 3}};

const slibreal_t Test_BremsstrahlungScreenedEmission::TEST_Z0[NTESTVALUES][6] = {
{1}, {2}, {1, 1}, {2, 1}, {3, 1, 1, 4}, {0, 3, 0, 1}, {0, 1, 1, 2, 1}, {1, 4, 4, 3, 1}, {3, 2, 1, 1, 0, 2}, {1, 0, 4, 2, 2, 2}};

const slibreal_t Test_BremsstrahlungScreenedEmission::TEST_DENS[NTESTVALUES][6] = {
{1.073500e+19}, 
{7.800400e+19}, 
{8.059500e+19, 7.451300e+19},
{6.523100e+19, 3.181100e+19},
{8.723700e+19, 3.280300e+19, 4.708600e+19, 9.909300e+19},
{9.563800e+19, 5.024300e+19, 4.180200e+19, 7.086400e+19},
{1.010200e+20, 1.009100e+20, 6.915900e+19, 4.325700e+19, 9.530600e+19},
{5.424000e+19, 1.004360e+20, 1.331800e+19, 6.324300e+19, 8.165000e+19},
{2.793000e+19, 4.365300e+19, 2.877100e+19, 4.219300e+19, 5.038600e+19, 6.485700e+19},
{1.487400e+19, 6.527300e+19, 3.748100e+19, 3.415000e+19, 3.431500e+19, 2.541600e+19}
};

struct TESTDATA {unsigned int nspecies; slibreal_t val;};
const struct TESTDATA TEST_NR_POW[10] = {
{1, 6.73045973479381e-10},
{1, 2.01662976950185e-09},
{2, 1.13293399609858e-09},
{2, 2.50881731824007e-09},
{4, 7.60986501112254e-09},
{4, 5.86883166300790e-09},
{5, 1.57595925024440e-08},
{5, 8.02396169262980e-09},
{6, 4.54479515488073e-09},
{6, 3.45373484297087e-09}
};

const slibreal_t Test_BremsstrahlungScreenedEmission::TEST_GAMMA_INTSPEC[NTESTVALUES][2] = {
{1,     0},
{5,     2.15962917566497e-09},
{10,	5.15999930299167e-09},
{15,	7.54998734239249e-09},
{25,	1.10284516497038e-08},
{30,	1.23582564705472e-08},
{40,	1.45255528103682e-08},
{50,	1.62554212133256e-08},
{75,	1.87447632097682e-08},
{100,	2.02556589173730e-08}
};

/**
 * Check that the formula for total emission
 * of synchrotron radiation has been implemented
 * correctly. Tests vs values obtained from same formulas implemented in matlab
 */
bool Test_BremsstrahlungScreenedEmission::CheckTotalEmission(const slibreal_t tol) {
    unsigned int i;
    Detector *det = GetDetector(0);
    slibreal_t pwr, corr, Delta;

    for (i = 0; i < NTESTVALUES; i++) {
	    slibreal_t Input_Z[6] = {TEST_Z[i][0], TEST_Z[i][1], TEST_Z[i][2], TEST_Z[i][3], TEST_Z[i][4], TEST_Z[i][5]};
        slibreal_t Input_DENS[6] = {TEST_DENS[i][0], TEST_DENS[i][1], TEST_DENS[i][2], TEST_DENS[i][3], TEST_DENS[i][4], TEST_DENS[i][5]};
        slibreal_t Input_Z0[6] = {TEST_Z0[i][0], TEST_Z0[i][1], TEST_Z0[i][2], TEST_Z0[i][3], TEST_Z0[i][4], TEST_Z0[i][5]};
        ConeBremsstrahlungScreenedEmission cbse(det, nullptr, TEST_NR_POW[i].nspecies, Input_Z, Input_Z0, Input_DENS);

        cbse.CalculateTotalEmission();
        pwr = cbse.GetTotalEmission();
        corr = TEST_NR_POW[i].val;
        
        Delta = fabs((pwr-corr)/corr);

        if (Delta >= tol) {
            this->PrintError("Total bremsstrahlung screened emission was not calculated correctly. Delta = %e", Delta);
            return false;
        }
    }
    return true;
}

/* 
 * Test so that the spectrum integrates to the correct value
 */
bool Test_BremsstrahlungScreenedEmission::CheckSpectrumEmission(const slibreal_t tol) {
    
    unsigned int i;
    Detector *det = GetDetector(50, 1, 50);
    MagneticFieldAnalytical2D *dummy_mf = GetMagneticField();
    slibreal_t pwr, corr, Delta;

    const unsigned int nspecies = TEST_NR_POW[NTESTVALUES-1].nspecies; 

    slibreal_t *Z = new slibreal_t[nspecies],
    *Z0 = new slibreal_t[nspecies], 
    *DENS = new slibreal_t[nspecies];
    for (i = 0; i < nspecies; i++){
        Z[i] = TEST_Z[NTESTVALUES-1][i];
        Z0[i] = TEST_Z0[NTESTVALUES-1][i];
        DENS[i] = TEST_DENS[NTESTVALUES-1][i];
    }
 
    for (i = 0; i < NTESTVALUES; i++) {
        ConeBremsstrahlungScreenedEmission cbse(det, nullptr, nspecies, Z, Z0, DENS);
        
        RadiationParticle *rp = GetRadiationParticle(i, det, dummy_mf);
        
        cbse.CalculateSpectrum(rp);
        cbse.IntegrateSpectrum();
        pwr = cbse.GetTotalEmission();
        corr = TEST_GAMMA_INTSPEC[i][IINTSPEC];
        
        Delta = fabs((pwr-corr)/corr);

        if (Delta >= tol) {
            this->PrintError("Total bremsstrahlung screened emission was not calculated correctly. Delta = %e", Delta);
            return false;
        }
    }
    return true;
}

/**
 * Returns a sample detector.
 */
Detector *Test_BremsstrahlungScreenedEmission::GetDetector(unsigned int nwavelengths, slibreal_t l0, slibreal_t l1) {
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
RadiationParticle *Test_BremsstrahlungScreenedEmission::GetRadiationParticle(unsigned int i, Detector *det, MagneticFieldAnalytical2D *mf) {
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
        gamma = TEST_GAMMA_INTSPEC[i][IGAMMA],
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
bool Test_BremsstrahlungScreenedEmission::Run(bool) {
    bool success = true;
    slibreal_t TOTEM_TOL=SQRT_REAL_EPSILON;
    slibreal_t SPECT_TOL=1e-3;

    if (CheckTotalEmission(TOTEM_TOL))
        this->PrintOK("Total emission is implemented correctly.");
    else
        success = false;

    if (CheckSpectrumEmission(SPECT_TOL))
        this->PrintOK("Spectrum integrates correctly.");
    else
        success = false;
    
    return success;
}

