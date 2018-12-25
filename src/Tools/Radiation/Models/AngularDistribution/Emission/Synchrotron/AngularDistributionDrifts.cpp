/**
 * Implementation of the angular distribution of
 * synchrotron radiation in first-order guiding-center
 * theory (including drift effects).
 */

#include <softlib/Integration/PeakedIntegration.h>
#include "Tools/Radiation/Models/AngularDistribution/ADSynchrotronEmission.h"

using namespace std;
using namespace __Radiation;

/**
 * FIRST (1) ORDER
 * Calculate the angular distribution of synchrotron
 * radiation in first-order guiding-center theory.
 */
slibreal_t ADSynchrotronEmission::CalculateAngularDistribution_FirstOrderGC(
    Vector<3> &n, slibreal_t __UNUSED__(sinMu), slibreal_t cosMu
) {
    gsl_function f;
    f.function = ADSynchrotronEmission::synchrotron_first_order_gc;
    f.params   = (void*)&angdistParams;

    // Observation vector
    angdistParams.n = -n;
    double abserr, pwr;

    gsl_integration_qags(
        &f, 0.0, 2.0*M_PI,
        this->qagsEpsAbs, this->qagsEpsRel,
        this->qagsLimit, this->qagsWorkspace,
        &pwr, &abserr
    );

    this->power = this->prefactor * pwr * (1.0 - this->beta*this->cosThetap*cosMu);
    return this->power;
}

/**
 * Internal routine for calculating the angular distribution
 * of synchrotron radiation in first-order guiding-center theory.
 * This routine evaluates the formula at a particular gyro-angle.
 */
double ADSynchrotronEmission::synchrotron_first_order_gc(double zeta, void *params) {
    struct angdist_params *p = (struct angdist_params*)params;
    Vector<3> betaDot, n_beta, n_betaXbetaDot, num, beta;
    double den, den2, den5;
    double cosZeta = cos(zeta),
        sinZeta = sqrt(1 - cosZeta*cosZeta);

    if (zeta > M_PI)
        sinZeta = -sinZeta;

    CalculateBeta(p, sinZeta, cosZeta, beta);
    CalculateBetaDot(p, sinZeta, cosZeta, beta, betaDot);

    _Rotate(beta, p->sinphi, p->cosphi);
    _Rotate(betaDot, p->sinphi, p->cosphi);

    n_beta = p->n - beta;

    Vector<3>::Cross(n_beta, betaDot, n_betaXbetaDot);
    Vector<3>::Cross(p->n, n_betaXbetaDot, num);

    den = 1.0 - beta.Dot(p->n);
    den2 = den * den;
    den5 = den2*den2*den;

    slibreal_t pwr = num.Dot(num)*p->betaDotFactor / den5;

    return pwr;
}

/**
 * Calculate the normalized particle velocity 'beta'
 * in first-order guiding-center theory.
 *
 * p:    Parameters to use when calculating 'beta'.
 * beta: Contains the normalized velocity on return.
 */
Vector<3> &ADSynchrotronEmission::CalculateBeta(
    struct angdist_params *p, double sinZeta,
    double cosZeta, Vector<3> &beta
) {
    slibreal_t
        sin2 = 2.0*sinZeta*cosZeta,
        cos2 = cosZeta*cosZeta - sinZeta*sinZeta;

    // First, calculate p
    beta = p->PGC + (sinZeta*p->sin + cosZeta*p->cos + sin2*p->sin2 + cos2*p->cos2);
    p->iGamma = 1.0 / sqrt(1.0 + beta.Dot(beta));
    beta *= p->iGamma;

    return beta;
}

/**
 * Calculate the particle's normalized acceleration,
 * 'beta-dot', in first-order guiding-center theory.
 * 
 * p:       Parameters to use when calculating the acceleration.
 * sin:     Sine of the gyro angle.
 * cos:     Cosine of the gyro angle.
 * beta:    The particle velocity, normalized to c.
 * betadot: Contains the acceleration on return.
 */
Vector<3> &ADSynchrotronEmission::CalculateBetaDot(
    struct angdist_params *p, slibreal_t sinZeta, slibreal_t cosZeta,
    const Vector<3> &beta, Vector<3> &betaDot
) {
    Vector<3> perpHat = -p->oneHat * sinZeta - p->twoHat * cosZeta;
    Vector<3> gradBrho= cosZeta*p->oneDotNablaB - sinZeta*p->twoDotNablaB;
    Vector<3>::Cross(beta, p->Bvec, betaDot);

    // Add 0'th order p times first-order B
    betaDot = betaDot + p->iGamma * p->rho * Vector<3>::Cross(p->ppar*p->bHat + p->pperp*perpHat, cosZeta*p->oneDotNablaB - sinZeta*p->twoDotNablaB);
    
    return betaDot;
}

/**
 * Initialize the current toroidal step.
 *
 * sinphi, cosphi: Sine/cosine of current toroidal angle.
 */
void ADSynchrotronEmission::InitializeToroidalStepAD(const slibreal_t sinphi, const slibreal_t cosphi) {
    if (!includeDrifts) return;

    angdistParams.sinphi = sinphi;
    angdistParams.cosphi = cosphi;
}

/**
 * Calculate quantities that do not depend on the
 * gyro angle.
 *
 * rp: Object representing particle emitting state.
 */
void ADSynchrotronEmission::PrepareFirstOrder(
    RadiationParticle *rp, struct angdist_params *p,
    slibreal_t *prefactor
) {
    const slibreal_t
        c = LIGHTSPEED;

    slibreal_t
        ppar = rp->GetPpar(),
        pperp = rp->GetPperp(),
        m = rp->GetMass(),
        q = rp->GetCharge(),
        B = rp->GetB(),
        **J = rp->GetBJacobian();

    p->B = B;
    p->m = m;
    p->q = q;
    p->gamma = rp->GetGamma();
    p->ppar = ppar;
    p->pperp = pperp;
    p->iGamma = 1.0 / p->gamma;
    p->gammapar2 = p->gamma*p->gamma / (1 + pperp*pperp);
    p->betaDotFactor = q*q / (p->gamma*p->gamma*m*m);

    p->sinphi = 0.0;
    p->cosphi = 1.0;

    // Guiding-center velocity (normalized to c)
    p->PGC = rp->GetP();

    *prefactor = q*q * p->gammapar2 / (16.0*M_PI*M_PI*EPS0*c);

    Vector<3> &bhat = p->bHat = rp->GetBhat();
    p->Bvec = rp->GetBvec();

    if (fabs(bhat[2]) > fabs(bhat[1]) && fabs(bhat[2]) > fabs(bhat[0])) {
        p->oneHat[0] = 1.0;
        p->oneHat[1] = 1.0;
        p->oneHat[2] = -(bhat[0] + bhat[1]) / bhat[2];
    } else if (fabs(bhat[1]) > fabs(bhat[0])) {
        p->oneHat[0] = 1.0;
        p->oneHat[1] = -(bhat[2] + bhat[0]) / bhat[1];
        p->oneHat[2] = 1.0;
    } else {
        p->oneHat[0] = -(bhat[1] + bhat[2]) / bhat[0];
        p->oneHat[1] = 1.0;
        p->oneHat[2] = 1.0;
    }

    p->oneHat.Normalize();
    Vector<3>::Cross(bhat, p->oneHat, p->twoHat);
    Vector<3> &one = p->oneHat, &two = p->twoHat;

    p->gradB = rp->GetGradB();
    p->curlB = rp->GetCurlB();
    p->curlBhat = rp->GetCurlBHat();

    slibreal_t rho = m*c*pperp / (fabs(q)*B);
    slibreal_t rhoPpar = rho*ppar;
    slibreal_t mMuQ = m*c*pperp*pperp / (2.0*q*B);  // In normalized units
    slibreal_t sgnq = (q>0?1.0:-1.0);

    p->rho = rho;

    Vector<3> kappa = Vector<3>::Cross(p->curlBhat, bhat);
    slibreal_t kappaDot1   = kappa.Dot(p->oneHat);
    slibreal_t kappaDot2   = kappa.Dot(p->twoHat);
    slibreal_t bGradLogB   = p->bHat.Dot(p->gradB) / B;
    slibreal_t oneGradLogB = p->oneHat.Dot(p->gradB) / B;
    slibreal_t twoGradLogB = p->twoHat.Dot(p->gradB) / B;

    // Gradients of individual components of B
    Vector<3> dBx(J[0]), dBy(J[1]), dBz(J[2]);
    // Corresponding gradients of individual components of bHat
    Vector<3> dbx = (dBx - bhat[0]*p->gradB) / B,
              dby = (dBy - bhat[1]*p->gradB) / B,
              dbz = (dBz - bhat[2]*p->gradB) / B;

    slibreal_t
        a1GradBhat =
            one[0]*one.Dot(dbx) +
            one[1]*one.Dot(dby) +
            one[2]*one.Dot(dbz) -
            two[0]*two.Dot(dbx) -
            two[1]*two.Dot(dby) -
            two[2]*two.Dot(dbz);
    slibreal_t
        a2GradBhat =
            one[0]*two.Dot(dbx) +
            one[1]*two.Dot(dby) +
            one[2]*two.Dot(dbz) +
            two[0]*one.Dot(dbx) +
            two[1]*one.Dot(dby) +
            two[2]*one.Dot(dbz);
    
    p->oneDotNablaB[0] = one.Dot(dBx);
    p->oneDotNablaB[1] = one.Dot(dBy);
    p->oneDotNablaB[2] = one.Dot(dBz);
    p->twoDotNablaB[0] = two.Dot(dBx);
    p->twoDotNablaB[1] = two.Dot(dBy);
    p->twoDotNablaB[2] = two.Dot(dBz);

    // Coefficients (by degree of sin(zeta)/cos(zeta))
    p->sin  = -pperp*sgnq*p->oneHat
            - rhoPpar*kappaDot2*p->bHat
            + 0.5*rhoPpar*bGradLogB*p->twoHat
            - 0.25*rhoPpar*(
                a2GradBhat*p->oneHat -
                a1GradBhat*p->twoHat
            );

    p->cos  = rhoPpar*kappaDot1*p->bHat
            - pperp*sgnq*p->twoHat
            - 0.5*rhoPpar*bGradLogB*p->oneHat
            + 0.25*rhoPpar*(
                a1GradBhat*p->oneHat +
                a2GradBhat*p->twoHat
            );

    p->sin2 = 0.5*mMuQ*a1GradBhat*p->bHat
            - 0.5*q*B*rho*rho/(m*c)*(
                oneGradLogB*p->oneHat -
                twoGradLogB*p->twoHat
            );
    p->cos2 = 0.5*mMuQ*a2GradBhat*p->bHat
            - 0.5*q*B*rho*rho/(m*c)*(
                twoGradLogB*p->oneHat +
                oneGradLogB*p->twoHat
            );
}

/**
 * Rotate the vector v in the xy-plane.
 *
 * v: Vector to rotate.
 * s: Sine of rotation angle.
 * c: Cosine of rotation angle.
 */
void ADSynchrotronEmission::_Rotate(Vector<3> &v, const slibreal_t s, const slibreal_t c) {
    slibreal_t a, b;
    a = v[0];
    b = v[1];

    v[0] = a*c + b*s;
    v[1] =-a*s + b*c;
}

/**
 * Rotate the gyro plane coordinate system 1-2. This
 * function is only intended for testing the implementation.
 *
 * rzeta: Angle by which to rotate the coordinate system.
 */
void ADSynchrotronEmission::RotateOneTwo(const slibreal_t rzeta) {
    struct angdist_params *p = &angdistParams;
    Vector<3> _1(p->oneHat), _2(p->twoHat);

    p->oneHat = _1*cos(rzeta) + _2*sin(rzeta);
    p->twoHat =-_1*sin(rzeta) + _2*cos(rzeta);
}

