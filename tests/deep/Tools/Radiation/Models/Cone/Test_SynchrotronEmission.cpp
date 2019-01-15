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
    ITHETAP = 1,
    IB = 2;

const unsigned int Test_SynchrotronEmission::NTESTPARTICLES = 15;
const slibreal_t Test_SynchrotronEmission::TESTPARTICLES[NTESTPARTICLES][3] = {
    // gamma,  thetap,  B
    { 2.0,     0.10,    5.0},
    { 2.0,     0.40,    5.0},
    { 2.0,     0.90,    5.0},

    { 5.0,     0.15,    5.0},
    { 5.0,     0.45,    5.0},
    { 5.0,     0.95,    5.0},

    {15.0,     0.20,    5.0},
    {15.0,     0.50,    5.0},
    {15.0,     1.00,    5.0},

    {55.0,     0.25,    5.0},
    {55.0,     0.55,    5.0},
    {55.0,     1.05,    5.0},

    {100.0,    0.30,    5.0},
    {100.0,    0.60,    5.0},
    {100.0,    1.10,    5.0}
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
    ConeSynchrotronEmission cse(det);
    slibreal_t pwr, corr, Delta;

    for (i = 0; i < NTESTPARTICLES; i++) {
        RadiationParticle *rp = GetRadiationParticle(i, det);
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
    slibreal_t pwr, corr, Delta, lambdac, gamma2, pperp2, gammapar2, betaperp2;
    Detector *dummyDet = GetDetector(0);

    for (i = 0; i < NTESTPARTICLES; i++) {
        RadiationParticle *rp = GetRadiationParticle(i, dummyDet);

        gamma2 = rp->GetGamma()*rp->GetGamma();
        pperp2 = rp->GetPperp()*rp->GetPperp();
        gammapar2 = gamma2 / (1.0 + pperp2);
        betaperp2 = pperp2 / gamma2;

        lambdac = GetLambdaC(rp);

        Detector *det = GetDetector(40001, lambdac/10.0, lambdac*1000.0);
        ConeSynchrotronEmission cse(det);

        cse.CalculateSpectrum(rp);
        cse.IntegrateSpectrum();
        // The following is necessary to get agreement at low energies
        pwr = cse.GetTotalEmission() * gammapar2 * betaperp2;
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
        l = 4.0*M_PI*ELECTRON_MASS*LIGHTSPEED / (3.0*ELEMENTARY_CHARGE),
        gamma = rp->GetGamma(),
        pperp = rp->GetPperp(),
        gammapar = gamma / sqrt(1.0 + pperp*pperp),
        B = rp->GetB();

    return (l*gammapar / (gamma*gamma*B));
}

/**
 * Returns the pre-defined RadiationParticle object
 * with index 'i'.
 *
 * i: Index in table 'predef_particles' with pre-defined
 *    particle parameters to test.
 */
RadiationParticle *Test_SynchrotronEmission::GetRadiationParticle(unsigned int i, Detector *det) {
    if (i >= NTESTPARTICLES)
        throw SOFTException("Trying to access non-existant test-particle.");

    slibreal_t
        Jdtdrho = 1.0,
        gamma = TESTPARTICLES[i][IGAMMA],
        p2 = gamma*gamma - 1.0,
        p  = sqrt(p2),
        cosThetap = cos(TESTPARTICLES[i][ITHETAP]),
        sinThetap = sqrt(1.0 - cosThetap*cosThetap),
        ppar = p * cosThetap,
        pperp = p * sinThetap,
        B = TESTPARTICLES[i][IB],
        m = ELECTRON_MASS, q = -ELEMENTARY_CHARGE,
        BB[3] = {1.0,0.0,0.0},
        Pp[3] = {p,0.0,0.0},
        Xx[3] = {0.0,0.0,0.0};
    Vector<3> P(Pp), X(Xx), Bvec(BB), bHat = Bvec / B;

    return new RadiationParticle(
        X, P,
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

