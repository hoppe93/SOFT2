/**
 * Implementation of the plain angular distribution
 * of synchrotron radiation (integrated over all wavelengths).
 */

#include <gsl/gsl_integration.h>
#include "Tools/Radiation/Models/AngularDistribution/ADSynchrotronEmission.h"

using namespace std;
using namespace __Radiation;

/**
 * Prepare to compute the angular distribution
 * of synchrotron radiation.
 *
 * rp: Particle emitting state.
 */
const slibreal_t ADSynchrotronEmission::ANGDIST_PREFAC_ND =
    1.0 / (8.0*M_PI*EPS0*LIGHTSPEED);
void ADSynchrotronEmission::PrepareAngularDistribution(RadiationParticle *rp) {
    slibreal_t
        q           = rp->GetCharge(),
        m           = rp->GetMass(),
        q4          = q*q*q*q,
        m2          = m*m;

    this->B         = rp->GetB();
    this->p         = sqrt(rp->GetP2());
    this->pperp     = rp->GetPperp();
    this->gamma     = rp->GetGamma();
    this->igamma    = 1.0 / gamma;
    this->igamma2   = igamma*igamma;
    this->beta      = p * igamma;
    this->sinThetap = pperp / p;
    this->cosThetap = fabs(rp->GetPpar() / p);

    if (!includeDrifts) {
        this->prefactor = ANGDIST_PREFAC_ND * q4 * B*B * beta*beta*sinThetap*sinThetap / (m2*(1.0+pperp*pperp));
    } else {
        ADSynchrotronEmission::PrepareFirstOrder(rp, &angdistParams, &(this->prefactor));
    }
}

/**
 * Calculate the plain angular distribution of synchrotron radiation.
 */
void ADSynchrotronEmission::CalculateAngularDistribution(
    Vector<3> &n, slibreal_t sinMu, slibreal_t cosMu
) {
    if (!includeDrifts)
        CalculateAngularDistribution_ZerothOrderGC(sinMu, cosMu);
    else
        CalculateAngularDistribution_FirstOrderGC(n, sinMu, cosMu);
}

/**
 * ZEROTH (0) ORDER
 * Calculate the angular distribution of synchrotron
 * radiation in zeroth-order guiding-center theory,
 * as given by Eq. (12) in [Hoppe et al., NF 2 (2018)].
 */
slibreal_t ADSynchrotronEmission::CalculateAngularDistribution_ZerothOrderGC(
    slibreal_t sinMu, slibreal_t cosMu
) {
    slibreal_t
        kappa = 1.0 - beta*cosMu*cosThetap,
        ikappa = 1.0 / kappa,
        ikappa2 = ikappa*ikappa,
        eta = ikappa * beta*sinMu*sinThetap,
        factor = prefactor * ikappa2;
    
    if (eta == 0.0) return factor;

    slibreal_t
        xi2 = 1.0/(1.0 - eta*eta),
        xi3 = xi2 * sqrt(xi2),
        xi5 = xi3*xi2,
        xi7 = xi5*xi2,

        p1  = 1.500 * xi5 - 0.500 * xi3,
        p2  = 0.625 * xi7 - 0.125 * xi5,
        p2f = sinMu*sinMu * igamma2 * ikappa2;

    this->power = factor * (p1 - p2f*p2);
    return this->power;
}

