/**
 * Implementation of the angular and spectral distribution
 * of synchrotron radiation, including polarization components.
 */

#include "Tools/Radiation/Models/AngularDistribution/ADSynchrotronEmission.h"
#include "Tools/Radiation/synchrotron_func.h"

using namespace std;
using namespace __Radiation;

/*************************
 * PREPARATION FUNCTIONS *
 *************************/
void ADSynchrotronEmission::PreparePolarization(RadiationParticle*) { }
void ADSynchrotronEmission::PrepareSpectrum(RadiationParticle *rp) {
    this->gamma = rp->GetGamma();
    this->igamma2 = 1.0 / (gamma*gamma);
    this->gamma3 = gamma*gamma*gamma;

    slibreal_t
        m = rp->GetMass(),
        q = rp->GetCharge(),
        B = rp->GetB(),
        p2 = rp->GetP2(),
        beta2 = p2 / (1.0 + p2),
        betapar2 = rp->GetPpar()*rp->GetPpar() / (1.0 + p2),
        gammapar = 1.0 / sqrt(1.0 - betapar2);

    this->prefactor = q*q*q*beta2*B / (16.0*M_PI*EPS0*gamma*m);
    this->beta = sqrt(beta2);

    this->lambdac = 4.0*M_PI*LIGHTSPEED*gammapar*igamma2*m / (3.0*q*B);
}

/**
 * Wrapper for calculating polarization of
 * synchrotron radiation.
 */
void ADSynchrotronEmission::CalculatePolarization(
    Vector<3> &n, slibreal_t sinMu,  slibreal_t cosMu
) { __CalculateSpectrum<true>(n, sinMu, cosMu); }

/**
 * Wrapper for calculating angular and spectral
 * distribution of synchrotron radiation.
 */
void ADSynchrotronEmission::CalculateSpectrum(
    Vector<3> &n, slibreal_t sinMu,  slibreal_t cosMu
) { __CalculateSpectrum<false>(n, sinMu, cosMu); }

/**
 * Function for calculating the spectrum of emitted synchrotron
 * radiation, including its polarization components.
 * 
 * calculatePolarization: If true, calculates the polarization
 *                        components of synchrotron radiation.
 * rp:                    Object representing the particle
 *                        emitting state.
 */
template<bool calculatePolarization>
void ADSynchrotronEmission::__CalculateSpectrum(
    Vector<3>&, slibreal_t sinMu,  slibreal_t cosMu
) {
    unsigned int i;
    slibreal_t
        cosPsi = cosMu*cosThetap + sinMu*sinThetap,
        sinPsi2 = 1.0 - cosPsi*cosPsi,
        mcospsi = 1.0-beta*cosPsi,
        xifac = gamma3*lambdac*sqrt(mcospsi*mcospsi*mcospsi/(0.5*beta*cosPsi));

    for (i = 0; i < nwavelengths; i++) {
        slibreal_t l = wavelengths[i];
        slibreal_t pfac = this->prefactor/(l*l * beta*cosPsi*(1.0 - beta*cosPsi));
        slibreal_t fac13 = (0.5*beta*cosPsi*sinPsi2)/(1.0 - beta*cosPsi);
        slibreal_t xi = xifac / l;

        slibreal_t
            xK13 = synchrotron_func3(xi),
            xK23 = synchrotron_func2(xi);

        I[i] = pfac * (xK23*xK23 + fac13*xK13*xK13);
    }
}

