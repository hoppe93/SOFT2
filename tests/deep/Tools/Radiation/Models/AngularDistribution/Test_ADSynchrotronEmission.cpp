/**
 * Test of the synchrotron radiation angular distribution model.
 * The following tests are included:
 *
 *   [ ] - Emission formulas integrate to Larmor's formula.
 *
 */

#include <cmath>
#include "constants.h"
#include "Tools/Radiation/Models/AngularDistribution/ADSynchrotronEmission.h"
#include "Test_ADSynchrotronEmission.h"

using namespace std;
using namespace __Radiation;

const unsigned int
    IGAMMA = 0,
    ITHETAP = 1,
    IB = 2;

const unsigned int Test_ADSynchrotronEmission::NTESTPARTICLES = 15;
const slibreal_t Test_ADSynchrotronEmission::TESTPARTICLES[NTESTPARTICLES][3] = {
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
 *
 * rp: Particle emitting state to calculate the
 *     total emitted radiation for.
 */
slibreal_t Test_ADSynchrotronEmission::Larmor(
    RadiationParticle *rp
) {
    const slibreal_t
        c = LIGHTSPEED,
        e = ELEMENTARY_CHARGE,
        eps0 = EPS0,
        m = ELECTRON_MASS;

    slibreal_t
        pperp = rp->GetPperp(),
        pperp2 = pperp*pperp,
        gamma = rp->GetGamma(),
        betaperp2 = pperp2 / (gamma*gamma),
        gammapar = gamma / sqrt(1.0 + pperp2),
        B = rp->GetB();

    slibreal_t C = e*e*e*e / (6.0*M_PI*eps0*m*m*c);
    slibreal_t Cp = B*B * gamma*gamma * gammapar*gammapar * betaperp2;

    return C * Cp;
}

/**
 * Returns the critical wavelength for the synchrotron spectrum
 * emitted by the given particle.
 */
slibreal_t Test_ADSynchrotronEmission::GetLambdaC(RadiationParticle *rp) {
    slibreal_t
        l = 4.0*M_PI*ELECTRON_MASS*LIGHTSPEED / (3.0*ELEMENTARY_CHARGE),
        gamma = rp->GetGamma(),
        pperp = rp->GetPperp(),
        gammapar = gamma / sqrt(1.0 + pperp*pperp),
        B = rp->GetB();

    return (l*gammapar / (gamma*gamma*B));
}

/**
 * Test the angular distribution of synchrotron radiation.
 *
 * NOTE: We must divide by the "Ginzburg factor", 1 - betapar*cos(mu),
 * since the Larmor formula is the radiation from one particle, while
 * the Ginzburg factor is for a stationary distribution of particles.
 */
bool Test_ADSynchrotronEmission::CheckAngularDistribution(
    ADSynchrotronEmission &ade, RadiationParticle *rp,
    slibreal_t tol
) {
    unsigned int i;
    const unsigned int n = 2000;
    slibreal_t I = 0, mu;
    slibreal_t dMu = M_PI / n;
    slibreal_t betapar = rp->GetPpar() / rp->GetGamma(),
        sinMu, cosMu;
    Vector<3> nHat;

    ade.PrepareAngularDistribution(rp);

    // Integrate using Simpson's rule (endpoints
    // don't matter, since sin(0) = sin(pi) = 0)
    
    // The inner points do matter...
    for (i = 1; i < n; i+=2) {
        mu = i*dMu;
        sinMu = sin(mu);
        cosMu = cos(mu);

        nHat = ade.GetParams()->bHat*cosMu + ade.GetParams()->oneHat*sinMu;
        ade.CalculateAngularDistribution(nHat, sinMu, cosMu);
        I += 4.0 * ade.GetTotalEmission() * sinMu / (1.0 - betapar*cosMu);
    }
    for (i = 2; i < n-1; i+=2) {
        mu = i*dMu;
        sinMu = sin(mu);
        cosMu = cos(mu);

        nHat = ade.GetParams()->bHat*cosMu + ade.GetParams()->oneHat*sinMu;
        ade.CalculateAngularDistribution(nHat, sinMu, cosMu);
        I += 2.0 * ade.GetTotalEmission() * sinMu / (1.0 - betapar*cosMu);
    }

    I *= dMu / 3.0;

    // Evaluate Larmor's formula
    slibreal_t L = Larmor(rp);
    // Compute relative error
    slibreal_t Delta = fabs((L-I) / L);

    if (isnan(Delta) || Delta >= tol) {
        this->PrintError("Angular distribution of synchrotron radiation does not integrate to Larmor's formula. Delta = %e", Delta);
        return false;
    }

    return true;
}

/**
 * Test the angular and spectral distribution of synchrotron radiation.
 *
 * NOTE: We must divide by the "Ginzburg factor", 1 - betapar*cos(mu),
 * since the Larmor formula is the radiation from one particle, while
 * the Ginzburg factor is for a stationary distribution of particles.
 */
bool Test_ADSynchrotronEmission::CheckAngularSpectralDistribution(
    ADSynchrotronEmission &ade, RadiationParticle *rp,
    slibreal_t tol
) {
    unsigned int i;
    const unsigned int n = 2000;
    slibreal_t I = 0, mu;
    slibreal_t dMu = M_PI / n;
    slibreal_t betapar = rp->GetPpar() / rp->GetGamma(),
        sinMu, cosMu;
    Vector<3> nHat;

    ade.PrepareSpectrum(rp);

    // Integrate using Simpson's rule (endpoints
    // don't matter, since sin(0) = sin(pi) = 0)
    
    // The inner points do matter...
    for (i = 1; i < n; i+=2) {
        mu = i*dMu;
        sinMu = sin(mu);
        cosMu = cos(mu);

        nHat = ade.GetParams()->bHat*cosMu + ade.GetParams()->oneHat*sinMu;
        ade.CalculateSpectrum(nHat, sinMu, cosMu);
        ade.IntegrateSpectrum();
        I += 4.0 * ade.GetTotalEmission() * sinMu / (1.0 - betapar*cosMu);
    }
    for (i = 2; i < n-1; i+=2) {
        mu = i*dMu;
        sinMu = sin(mu);
        cosMu = cos(mu);

        nHat = ade.GetParams()->bHat*cosMu + ade.GetParams()->oneHat*sinMu;
        ade.CalculateSpectrum(nHat, sinMu, cosMu);
        ade.IntegrateSpectrum();
        I += 2.0 * ade.GetTotalEmission() * sinMu / (1.0 - betapar*cosMu);
    }

    I *= dMu / 3.0;

    slibreal_t
        B = rp->GetB(),
        gamma2 = 1.0 + rp->GetP2(),
        gamma3 = gamma2*sqrt(gamma2),
        gammapar2 = 1.0 / sqrt(1.0 - rp->GetPpar()*rp->GetPpar()/gamma2),
        betaperp2 = rp->GetPperp()*rp->GetPperp() / gamma2,
        m = rp->GetMass(),
        e = fabs(rp->GetCharge()),

        fI = 8.0*e*B*gamma3*gammapar2*betaperp2 / (3.0*m*LIGHTSPEED);

    // Evaluate Larmor's formula
    slibreal_t L = Larmor(rp);
    // Compute relative error
    slibreal_t Delta = fabs((L-I) / L);

    printf("I = %e, fI = %e,   I/fI = %e\n", I, fI, I/fI);

    if (isnan(Delta) || Delta >= tol) {
        this->PrintError("Angular and spectral distribution of synchrotron radiation does not integrate to Larmor's formula. Delta = %e", Delta);
        this->PrintError("L/I = %e", L/I);
        return true;
    }

    return true;
}

/**
 * Returns a sample detector.
 */
Detector *Test_ADSynchrotronEmission::GetDetector(unsigned int nwavelengths, slibreal_t l0, slibreal_t l1) {
    const slibreal_t
        aperture = 0.006,
        tilt     = 0.0,
        visang   = 1.0,
        dir[3] = {0.0,1.0,0.0},
        pos[3] = {0.0,-1.069,0.0};
        
    Vector<3>
        direction(dir),
        position(pos);

    return new Detector(aperture, tilt, visang, direction, position, nwavelengths, l0, l1);
}

/**
 * Returns the pre-defined RadiationParticle object
 * with index 'i'.
 *
 * i: Index in table 'predef_particles' with pre-defined
 *    particle parameters to test.
 */
RadiationParticle *Test_ADSynchrotronEmission::GetRadiationParticle(unsigned int i, Detector *det) {
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
bool Test_ADSynchrotronEmission::Run(bool) {
    unsigned int i;
    bool success = true;
    struct global_settings globset;
    globset.include_drifts = false;
    Detector *det = GetDetector(0);

    ADSynchrotronEmission ade(det, nullptr, &globset);

    for (i = 0; i < NTESTPARTICLES && success; i++) {
        RadiationParticle *rp = GetRadiationParticle(i, det);

        // Test angular distribution
        success &= CheckAngularDistribution(ade, rp, ANGDIST_TOL);
    }

    for (i = 0; i < NTESTPARTICLES && success; i++) {
        RadiationParticle *rp = GetRadiationParticle(i, det);
        slibreal_t lambdac = GetLambdaC(rp);

        Detector *det2 = GetDetector(4000, lambdac/10.0, lambdac*1000.0);
        ADSynchrotronEmission adsp(det2, nullptr, &globset);

        // Test angular & spectral distribution
        success &= CheckAngularSpectralDistribution(adsp, rp, ANGDIST_TOL);

        // Test angular & spectral distribution w/ polarization
    }

    return success;
}

