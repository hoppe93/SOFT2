/**
 * Implementation of the angular and spectral distribution
 * of synchrotron radiation, including polarization components.
 */

#include "Tools/Radiation/Models/AngularDistribution/ADSynchrotronEmission.h"

using namespace std;
using namespace __Radiation;

/*************************
 * PREPARATION FUNCTIONS *
 *************************/
void ADSynchrotronEmission::PreparePolarization(RadiationParticle *rp) { }
void ADSynchrotronEmission::PrepareSpectrum(RadiationParticle *rp) { }

/**
 * Wrapper for calculating polarization of
 * synchrotron radiation.
 */
slibreal_t ADSynchrotronEmission::CalculatePolarization(
    Vector<3> &n, slibreal_t sinMu,  slibreal_t cosMu
) { return __CalculateSpectrum<true>(n, sinMu, cosMu); }

/**
 * Wrapper for calculating angular and spectral
 * distribution of synchrotron radiation.
 */
slibreal_t ADSynchrotronEmission::CalculateSpectrum(
    Vector<3> &n, slibreal_t sinMu,  slibreal_t cosMu
) { return __CalculateSpectrum<false>(n, sinMu, cosMu); }

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
slibreal_t ADSynchrotronEmission::__CalculateSpectrum(
    Vector<3> &n, slibreal_t sinMu,  slibreal_t cosMu
) {
    throw ADSynchrotronException("Not implemented yet.");
}

