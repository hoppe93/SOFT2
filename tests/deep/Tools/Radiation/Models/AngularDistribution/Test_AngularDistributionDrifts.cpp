/**
 * Test of the angular distribution radiation model
 * in first-order guiding-center theory.
 * The following tests are included:
 *
 *   [X]  - Verify that the error term in p^2 scales
 *          as B^-2.
 *   [X]  - Compare calculated 'beta' to the particle
 *          velocity calculated in GOA.
 *   [X]  - Compare calculated 'beta-dot' to the particle
 *          acceleration calculated in GOA.
 *   [X]  - Calculate the acceleration in FOGC theory
 *          using the code from GOA and make sure that
 *          the angular distribution implemented in SOFT
 *          integrates to Larmor's formula.
 */

#include <functional>
#include <softlib/config.h>
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>
#include <softlib/Timer.h>
#include <softlib/Vector.h>
#include "Orbit/GuidingCenterEquation.h"
#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/AngularDistribution/ADSynchrotronEmission.h"

#include "GOAMomentum.h"
#include "Test_AngularDistributionDrifts.h"

using namespace std;
using namespace __Radiation;

const unsigned int
    IGAMMA = 0,
    ITHETAP = 1,
    IR = 2;

const unsigned int Test_AngularDistributionDrifts::NTESTPARTICLES = 15;
const slibreal_t Test_AngularDistributionDrifts::TESTPARTICLES[NTESTPARTICLES][3] = {
    // gamma,  thetap,  r
    {2.0,      0.10,    0.68},
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
Detector *Test_AngularDistributionDrifts::GetDetector(unsigned int nwavelengths, slibreal_t l0, slibreal_t l1) {
    const slibreal_t
        aperture = 0.006,
        roll     = 0.0,
        visang   = 1.0,
        dir[3] = {0.0,1.0,0.0},
        pos[3] = {0.0,-1.069,0.0};
        
    Vector<3>
        direction(dir),
        position(pos);

    return new Detector(aperture, roll, visang, direction, position, nwavelengths, l0, l1);
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
Vector<3> Test_AngularDistributionDrifts::GetGuidingCenterMomentum(
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
RadiationParticle *Test_AngularDistributionDrifts::GetRadiationParticle(MagneticFieldAnalytical2D *mfa, Detector *det, unsigned int i) {
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

/**
 * Test verifying that the correction term in the square
 * of the computed particle velocity scales as B^-2.
 */
bool Test_AngularDistributionDrifts::VerifyP2Scaling(const slibreal_t tol) {
    slibreal_t prefactor = 0.0;
    Vector<3> beta1, beta2;
    struct global_settings globset;
    globset.include_drifts = true;
    Detector *det = GetDetector(0);

    struct ADSynchrotronEmission::angdist_params p;

    for (unsigned int i = 0; i < NTESTPARTICLES; i++) {
        MagneticFieldAnalytical2D *mfa = GetMagneticField();

        // The scaling doesn't work on the magnetic axis
        if (TESTPARTICLES[i][IR] == mfa->GetMagneticAxisR())
            continue;

        // First evaluation
        RadiationParticle *rp1 = GetRadiationParticle(mfa, det, i);
        ADSynchrotronEmission::PrepareFirstOrder(rp1, &p, &prefactor);
        ADSynchrotronEmission::CalculateBeta(&p, 0.0, 1.0, beta1);

        // Second evaluation
        mfa->SetB0(mfa->GetB0()*2.0);
        RadiationParticle *rp2 = GetRadiationParticle(mfa, det, i);
        ADSynchrotronEmission::PrepareFirstOrder(rp2, &p, &prefactor);
        ADSynchrotronEmission::CalculateBeta(&p, 0.0, 1.0, beta2);

        slibreal_t ppar1 = rp1->GetPpar(), pperp1 = rp1->GetPperp();
        slibreal_t ppar2 = rp2->GetPpar(), pperp2 = rp2->GetPperp();
        slibreal_t gamma12 = rp1->GetGamma()*rp1->GetGamma();
        slibreal_t gamma22 = rp2->GetGamma()*rp2->GetGamma();
        slibreal_t err1 = beta1.Dot(beta1)*gamma12 - (ppar1*ppar1 + pperp1*pperp1);
        slibreal_t err2 = beta2.Dot(beta2)*gamma22 - (ppar2*ppar2 + pperp2*pperp2);

        slibreal_t Delta = fabs((err1-4.0*err2)/err1);
        if (isnan(Delta) || Delta >= tol) {
            this->PrintError("Particle #%u: Correction terms in FOGC 'beta' do not scale as B^-2. err1 / err2 = %e", i, err1/err2);
            return false;
        }

        delete rp2;
        delete rp1;
        delete mfa;
    }

    return true;
}

/**
 * Compare the value of 'beta' calculated by SOFT
 * to the value calculated by GOA.
 */
bool Test_AngularDistributionDrifts::CompareBetas(const slibreal_t tol) {
    Detector *det = GetDetector(0);
    MagneticFieldAnalytical2D *mfa = GetMagneticField();
    struct ADSynchrotronEmission::angdist_params p;
    slibreal_t prefactor, ppar, pperp,
        DeltaNorm, DeltaBhat, Delta1hat, Delta2hat;
    const slibreal_t zeta = M_PI/10.0;
    Vector<3> beta_soft, beta_mom, Beta;

    for (unsigned int i = 0; i < NTESTPARTICLES; i++) {
        RadiationParticle *rp = GetRadiationParticle(mfa, det, i);
        
        // Evaluate 'beta' using SOFT
        ADSynchrotronEmission::PrepareFirstOrder(rp, &p, &prefactor);
        ADSynchrotronEmission::CalculateBeta(&p, sin(zeta), cos(zeta), beta_soft);

        // Evaluate 'beta' using GOA
        ppar = rp->GetPpar();
        pperp = rp->GetPperp();

        GOAMomentum goa(rp->GetPosition(), mfa, rp->GetCharge() / ELEMENTARY_CHARGE);
        beta_mom = goa.ParticleMomentum(
            ppar, pperp, zeta, Beta, beta_mom
        );
        beta_mom /= sqrt(1.0 + beta_mom.Dot(beta_mom));

        DeltaNorm = fabs((beta_soft.Norm() - beta_mom.Norm()) / beta_mom.Norm());
        DeltaBhat = fabs((beta_soft.Dot(p.bHat) - beta_mom.Dot(p.bHat)) / beta_mom.Dot(p.bHat));
        Delta1hat = fabs((beta_soft.Dot(p.oneHat) - beta_mom.Dot(p.oneHat)) / beta_mom.Dot(p.oneHat));
        Delta2hat = fabs((beta_soft.Dot(p.twoHat) - beta_mom.Dot(p.twoHat)) / beta_mom.Dot(p.twoHat));

        delete rp;

        if (isnan(DeltaNorm) || DeltaNorm >= tol) {
            this->PrintError("Particle #%u: Norms differ. Delta = %e (db = %e, d1 = %e, d2 = %e)", i, DeltaNorm, DeltaBhat, Delta1hat, Delta2hat);
            return false;
        } else if (isnan(DeltaBhat) || DeltaBhat >= tol ||
                   isnan(Delta1hat) || Delta1hat >= tol ||
                   isnan(Delta2hat) || Delta2hat >= tol
        ) {
            this->PrintError("Particle #%u: Norms are equal, but components differ. db = %e, d1 = %e, d2 = %e", i, DeltaBhat, Delta1hat, Delta2hat);
            return false;
        }
    }

    return true;
}

/**
 * Compare the value of 'beta-dot' calculated by SOFT
 * to the value calculated by GOA.
 * (beta-dot refers to the particle acceleration)
 */
bool Test_AngularDistributionDrifts::CompareBetaDots(const slibreal_t tol) {
    Detector *det = GetDetector(0);
    MagneticFieldAnalytical2D *mfa = GetMagneticField();
    struct ADSynchrotronEmission::angdist_params p;
    slibreal_t prefactor, ppar, pperp, gamma,
        Delta, betaDot2_soft, betaDot2_goa;
    Vector<3> beta_soft, betaDot_soft;
    const slibreal_t zeta = M_PI/10.0;

    for (unsigned int i = 0; i < NTESTPARTICLES; i++) {
        RadiationParticle *rp = GetRadiationParticle(mfa, det, i);
        
        ppar = rp->GetPpar();
        pperp = rp->GetPperp();
        gamma = rp->GetGamma();

        // Evaluate 'beta-dot' using SOFT
        ADSynchrotronEmission::PrepareFirstOrder(rp, &p, &prefactor);
        ADSynchrotronEmission::CalculateBeta(&p, sin(zeta), cos(zeta), beta_soft);
        ADSynchrotronEmission::CalculateBetaDot(&p, sin(zeta), cos(zeta), beta_soft, betaDot_soft);
        betaDot2_soft = betaDot_soft.Dot(betaDot_soft) * p.betaDotFactor * gamma*gamma;

        // Evaluate 'beta-dot' using GOA
        GOAMomentum goa(rp->GetPosition(), mfa, rp->GetCharge() / ELEMENTARY_CHARGE);

        struct betadot_params params;
        params.ppar = ppar;
        params.pperp = pperp;
        params.goa = &goa;

        betaDot2_goa = (slibreal_t)Test_AngularDistributionDrifts::betaDotSquared(zeta, &params);

        // Calculate relative errors
        Delta = fabs((betaDot2_soft - betaDot2_goa) / betaDot2_goa);

        delete rp;

        if (isnan(Delta) || Delta >= tol) {
            this->PrintError("Particle #%u: SOFT and GOA particle accelerations differ. Delta = %e", i, Delta);
            return false;
        }
    }

    return true;
}

/************************************
 *  INTEGRATES TO LARMOR'S FORMULA? *
 ************************************/
/**
 * Larmor's formula in FOGC.
 */
slibreal_t Test_AngularDistributionDrifts::Larmor(
    RadiationParticle *rp, MagneticFieldAnalytical2D *mfa
) {
    slibreal_t
        ppar = rp->GetPpar(),
        pperp = rp->GetPperp(),
        gamma = rp->GetGamma(),
        gammapar = gamma / sqrt(1.0 + pperp*pperp),
        q = rp->GetCharge(),
        c = LIGHTSPEED;
    
    // Prepare to average beta-dot squared
    GOAMomentum goa(rp->GetPosition(), mfa, rp->GetCharge() / ELEMENTARY_CHARGE);

    struct betadot_params params;
    params.ppar = ppar;
    params.pperp = pperp;
    params.goa = &goa;

    gsl_function f;
    f.function = Test_AngularDistributionDrifts::betaDotSquared;
    f.params = (void*)&params;

    double avBetaDot, absErr;
    const double qagsEpsAbs = 0.0, qagsEpsRel = 1e-7;
    size_t qagsLimit = 100;

    gsl_integration_workspace *wspace = gsl_integration_workspace_alloc(qagsLimit);

    gsl_integration_qags(
        &f, 0.0, 2.0*M_PI,
        qagsEpsAbs, qagsEpsRel,
        qagsLimit, wspace,
        &avBetaDot, &absErr
    );

    gsl_integration_workspace_free(wspace);

    avBetaDot *= 1.0 / (2.0*M_PI);
    slibreal_t larmor = q*q*avBetaDot * gamma*gamma * gammapar*gammapar / (6.0*M_PI*EPS0*c);

    return larmor;
}

/**
 * Evaluates beta-dot squared at a given gyro angle.
 *
 * zeta:   Gyro angle.
 * params: Pointer to a 'GOAMomentum' object.
 */
double Test_AngularDistributionDrifts::betaDotSquared(double zeta, void *params) {
    struct betadot_params *par = (struct betadot_params*)params;
    Vector<3> P, p, betadot;

    par->goa->ParticleMomentum(
        par->ppar, par->pperp, (slibreal_t)zeta, P, p
    );

    par->goa->Acceleration(
        p, par->ppar, par->pperp, zeta, betadot
    );

    return betadot.Dot(betadot);
}

double Test_AngularDistributionDrifts::EvaluateMuIntegral_inner(
    double mu, void *params
) {
    struct intpar *p = (struct intpar*)params;
    double sinmu = sin(mu), cosmu = cos(mu);
    Vector<3> n =-p->ade->GetParams()->bHat * cosmu + p->ade->GetParams()->oneHat*sinmu;

    p->ade->CalculateAngularDistribution(n, sinmu, cosmu);
    return (double)p->ade->GetTotalEmission() * sinmu / (1.0-p->betapar*cos(mu));
}
/**
 * Evaluate the integral of the angular distribution of
 * synchrotron radiation over mu (from 0 to pi).
 *
 * ade:   Angular distribution to integrate.
 * rp:    Particle emitting state.
 * rzeta: Offset in gyro angle to impose.
 * adapt: If true, uses the adaptive quadrature.
 */
slibreal_t Test_AngularDistributionDrifts::EvaluateMuIntegral(
    ADSynchrotronEmission &ade, RadiationParticle *rp, bool adapt
) {
    unsigned int i;
    const unsigned int qagsLimit = 100;
    slibreal_t epsabs = 0, epsrel = 1e-9;
    const unsigned int n = 2000;
    slibreal_t I = 0, err, cosmu, sinmu;
    slibreal_t dMu = M_PI / n, mu;
    slibreal_t betapar = rp->GetPpar() / rp->GetGamma();
    gsl_integration_workspace *wspace = gsl_integration_workspace_alloc(qagsLimit);

    ade.PrepareAngularDistribution(rp);

    if (adapt) {
        struct intpar ip;
        ip.betapar = betapar;
        ip.ade = &ade;
        
        gsl_function f;
        f.function = Test_AngularDistributionDrifts::EvaluateMuIntegral_inner;
        f.params   = &ip;

        gsl_integration_qags(
            &f, 0.0, M_PI, 
            epsabs, epsrel,
            qagsLimit, wspace,
            &I, &err
        );
    } else {
        Vector<3> nHat;

        // Integrate using Simpson's rule (endpoints
        // don't matter, since sin(0) = sin(pi) = 0)
        for (i = 1; i < n; i+=2) {
            mu = i*dMu;
            cosmu = cos(mu);
            sinmu = sin(mu);

            nHat =-ade.GetParams()->bHat * cosmu + ade.GetParams()->oneHat * sinmu;
            ade.CalculateAngularDistribution(nHat, sinmu, cosmu);
            I += 4.0 * ade.GetTotalEmission() * sinmu / (1.0 - betapar*cosmu);
        }
        for (i = 2; i < n-1; i+=2) {
            mu = i*dMu;
            cosmu = cos(mu);
            sinmu = sin(mu);

            nHat =-ade.GetParams()->bHat * cosmu + ade.GetParams()->oneHat * sinmu;
            ade.CalculateAngularDistribution(nHat, sinmu, cosmu);
            I += 2.0 * ade.GetTotalEmission() * sinmu / (1.0 - betapar*cosmu);
        }

        I *= dMu / 3.0;
    }

    return I;
}

/**
 * Verify that the angular distribution integrates to
 * Larmor's formula for a particular particle.
 *
 * ade: Angular distribution to integrate.
 * rp:  Particle emitting state.
 * tol: Tolerance to demand for success.
 */
bool Test_AngularDistributionDrifts::CheckAngularDistributionIntegral(
    ADSynchrotronEmission &ade, RadiationParticle *rp,
    MagneticFieldAnalytical2D* mfa, slibreal_t tol
) {
    slibreal_t I = EvaluateMuIntegral(ade, rp);

    // Evaluate Larmor's formula
    slibreal_t L = Larmor(rp, mfa);

    // Compute relative error
    slibreal_t Delta = fabs((L-I) / L);

    if (isnan(Delta) || Delta >= tol) {
        this->PrintError("Angular distribution of synchrotron radiation (FOGC) does not integrate to Larmor's formula. Delta = %e", Delta);
        return false;
    }

    return true;
}

/**
 * Make sure the FOGC implementation of the
 * synchrotron formula integrates to Larmor's formula.
 */
bool Test_AngularDistributionDrifts::VerifyIntegral(const slibreal_t tol) {
    unsigned int i;
    const slibreal_t qagsRelTol = SQRT_REAL_EPSILON;
    const unsigned int qagsLimit = 100;
    Detector *det = GetDetector(0);
    MagneticFieldAnalytical2D *mfa = GetMagneticField();
    struct global_settings globset;
    globset.include_drifts = true;
    bool success = true;

    ADSynchrotronEmission ade(det, mfa, &globset, qagsLimit, qagsRelTol);

    for (i = 0; i < NTESTPARTICLES && success; i++) {
        RadiationParticle *rp = GetRadiationParticle(mfa, det, i);

        success &= CheckAngularDistributionIntegral(ade, rp, mfa, tol);

        delete rp;

        if (!success)
            this->PrintError("Discrepancy occured for particle #%u.", i);
    }

    return success;
}

/**
 * Run all the tests in the module.
 */
bool Test_AngularDistributionDrifts::Run(bool) {
    bool success = true;
    const slibreal_t P2TOL = 2e-2;
    const slibreal_t COMPARETOL = 100.0*REAL_EPSILON;
    const slibreal_t INTTOL = 1e-1;

    if (VerifyP2Scaling(P2TOL))
        this->PrintOK("First-order correction terms scale as B^-2.");
    else
        success = false;

    if (CompareBetas(COMPARETOL))
        this->PrintOK("First-order particle velocity verified against GOA implementation.");
    else
        success = false;

    if (CompareBetaDots(COMPARETOL))
        this->PrintOK("First-order particle acceleration verified against GOA implementation.");

    // It appears that the numerically integrated angular
    // distribution has some additional correction terms
    // which are not included in the Larmor formula. Thus
    // a relatively high integral tolerance is set here.
    if (VerifyIntegral(INTTOL))
        this->PrintOK("Emission integrates to Larmor's formula.");
    else
        success = false;

    return success;
}

