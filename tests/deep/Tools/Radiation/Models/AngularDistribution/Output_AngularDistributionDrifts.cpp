/**
 * Output the angular distribution with drifts.
 */

#include <iostream>
#include <softlib/config.h>
#include <softlib/Vector.h>
#include "Orbit/GuidingCenterEquation.h"
#include "SOFT.h"
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/AngularDistribution/ADSynchrotronEmission.h"

#include "Output_AngularDistributionDrifts.h"

using namespace std;
using namespace __Radiation;

const unsigned int
    IGAMMA = 0,
    ITHETAP = 1,
    IR = 2;

const unsigned int Output_AngularDistributionDrifts::NTESTPARTICLES = 15;
const slibreal_t Output_AngularDistributionDrifts::TESTPARTICLES[NTESTPARTICLES][3] = {
    // gamma,  thetap,  r
    //{50.009999000199947,      0.2,    0.72},
    {25.0,      0.2,    0.72},
    {2.0,      0.30,    0.58},
    {2.0,      0.50,    0.78},
    {2.0,      0.70,    0.62},
    {2.0,      0.90,    0.72},

    {8.0,      0.10,    0.68},
    {8.0,      0.30,    0.58},
    {8.0,      0.50,    0.78},
    {8.0,      0.70,    0.62},
    {8.0,      0.90,    0.72},

    {20.0,     0.10,    0.68},
    {20.0,     0.30,    0.58},
    {20.0,     0.50,    0.78},
    {20.0,     0.70,    0.62},
    {20.0,     0.90,    0.72},
};

/**
 * Returns a sample detector.
 */
Detector *Output_AngularDistributionDrifts::GetDetector(unsigned int nwavelengths, slibreal_t l0, slibreal_t l1) {
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
 * Uses the 'GuidingCenterEquation' class to calculate the
 * guiding-center momentum vector.
 *
 * mf:    Pointer to magnetic field object used
 * mfd:   Magnetic field quantities evaluated in X
 * X:     Position of guiding-center
 * ppar:  Parallel momentum of particle
 * pperp: Perpendicular momentum of particle
 */
Vector<3> Output_AngularDistributionDrifts::GetGuidingCenterMomentum(
    MagneticField2D *mf, struct magnetic_field_data &mfd,
    const Vector<3> &X, const slibreal_t ppar, const slibreal_t pperp
) {
    struct global_settings globset;
    globset.include_drifts = true;
    GuidingCenterEquation gce(mf, &globset);
    Vector<6> dzdt, zval;
    slibreal_t 
        m = ELECTRON_MASS, q = -ELEMENTARY_CHARGE,
        c = LIGHTSPEED;

    zval[GuidingCenterEquation::COORD_X]    = X[0];
    zval[GuidingCenterEquation::COORD_Y]    = X[1];
    zval[GuidingCenterEquation::COORD_Z]    = X[2];
    zval[GuidingCenterEquation::COORD_PPAR] = ppar;
    zval[GuidingCenterEquation::COORD_MU]   = pperp*pperp*m*c*c / (2.0*mfd.Babs);
    zval[GuidingCenterEquation::COORD_ZETA] = 0.0;

    Particle p;
    p.SetCharge(q);
    p.SetMass(m);
    p.InitializeMomentum(
        Particle::COORDINATE_PPAR, Particle::COORDINATE_PPERP,
        ppar, pperp, 1.0, 1.0
    );
    p.InitializePosition(Particle::POSITION_GUIDINGCENTER, X[0], X[2], 1.0, 0.0);

    gce.InitializeParticle(&p, zval);
    gce.Evaluate(0.0, zval, dzdt);

    slibreal_t gamma = sqrt(1.0 + ppar*ppar + pperp*pperp);
    Vector<3> PP;
    PP[0] = dzdt[GuidingCenterEquation::COORD_X] * gamma / LIGHTSPEED;
    PP[1] = dzdt[GuidingCenterEquation::COORD_Y] * gamma / LIGHTSPEED;
    PP[2] = dzdt[GuidingCenterEquation::COORD_Z] * gamma / LIGHTSPEED;

    return PP;
}

/**
 * Constructs a sample radiation particle.
 *
 * mfa: Magnetic field in which the test particle should be created.
 * i:   Index of test particle to return.
 */
RadiationParticle *Output_AngularDistributionDrifts::GetRadiationParticle(MagneticFieldAnalytical2D *mfa, Detector *det, unsigned int i) {
    if (i >= NTESTPARTICLES)
        throw SOFTException("Trying to access non-existant test particle.");
    
    slibreal_t
        Jdtdrho = 1.0,
        gamma = TESTPARTICLES[i][IGAMMA],
        p2 = gamma*gamma - 1.0,
        p = sqrt(p2),
        cosThetap = cos(TESTPARTICLES[i][ITHETAP]),
        sinThetap = sqrt(1.0 - cosThetap*cosThetap),
        ppar = p*cosThetap,
        pperp = p*sinThetap,
        m = ELECTRON_MASS, q = -ELEMENTARY_CHARGE,
        XX[3] = {TESTPARTICLES[i][IR], 0.0, 0.0};

    MagneticField2D *mf = (MagneticField2D*)mfa;
    struct magnetic_field_data mfd = mf->EvalDerivatives(XX);

    // Use the Orbit module to evaluate the GC velocity
    Vector<3> PP = GetGuidingCenterMomentum(
        mf, mfd, XX, ppar, pperp
    );

    Vector<3> B(mfd.B), bHat(B/mfd.Babs);
    
    return new RadiationParticle(
        XX, PP, Jdtdrho, ppar, pperp,
        gamma, p2, det->GetPosition(),
        mfd.Babs, B, bHat, m, q, 0, 0, 0,
        mfd.gradB, mfd.curlB, mfd.J
    );
}

void Output_AngularDistributionDrifts::OutputIntegrand() {
    struct ADSynchrotronEmission::angdist_params p;
    struct global_settings globset;
    globset.include_drifts = true;
    slibreal_t prefactor, dZeta;
    unsigned int i;
    const unsigned int nParticle = 0,
        nZeta = 400;
    MagneticFieldAnalytical2D *mfa = GetMagneticField();
    Detector *det = GetDetector(0);

    dZeta = 2.0*M_PI / (nZeta-1);
    slibreal_t integrand[nZeta], zeta[nZeta];
    for (i = 0; i < nZeta; i++)
        zeta[i] = i*dZeta;

    RadiationParticle *rp = GetRadiationParticle(mfa, det, nParticle);
    slibreal_t
        gamma = rp->GetGamma(),
        //cosMu = rp->GetPpar() / sqrt(rp->GetP2()),
        cosMu = 0.98852765918194963,
        thetap = acos(cosMu),
        sinMu = sqrt(1.0 - cosMu*cosMu),
        cosDrift = rp->GetP().Dot(rp->GetBhat()) / rp->GetP().Norm();

    ADSynchrotronEmission::PrepareFirstOrder(rp, &p, &prefactor);
    p.n = cosMu * p.bHat + sinMu * p.oneHat;

    for (i = 0; i < nZeta; i++)
        integrand[i] = (slibreal_t)ADSynchrotronEmission::synchrotron_first_order_gc(zeta[i], &p);

    // Output
    for (i = 0; i < nZeta; i++)
        cerr << gamma << "," << thetap << "," << cosDrift << "," << zeta[i] << "," << integrand[i] << endl;
}

bool Output_AngularDistributionDrifts::Run(bool runningAll) {
    if (!runningAll)
        OutputIntegrand();

    return true;
}

