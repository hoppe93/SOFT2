/**
 * Module for evaluating synchrotron emission in
 * the cone model.
 */

#include "Tools/Radiation/Models/Cone/ConeUnitEmission.h"
#include "Tools/Radiation/RadiationParticle.h"

using namespace __Radiation;

/**
 * Evaluates the total emission and spectrum of
 * the unit emission.
 *
 * rp:           Object representing particle emitting state (unused).
 * polarization: If true, evaluate Stokes parameters (unused).
 */
void ConeUnitEmission::HandleParticle(RadiationParticle*, bool) {
    CalculateTotalEmission();
    if (nwavelengths > 0)
        CalculateSpectrum();
}

/**
 * Set total emission.
 */
void ConeUnitEmission::CalculateTotalEmission() {
    this->power = 1.0;
}

/**
 * Set spectral emission.
 */
void ConeUnitEmission::CalculateSpectrum() {
    for (unsigned int i = 0; i < nwavelengths; i++)
        I[i] = 1.0;
}

