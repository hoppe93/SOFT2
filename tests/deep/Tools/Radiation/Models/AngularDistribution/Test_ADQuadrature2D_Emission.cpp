/**
 * Implementation of the emission class for the
 * ADQuadrature2D test module.
 */

#include <softlib/config.h>
#include <softlib/Vector.h>
#include "Tools/Radiation/RadiationParticle.h"
#include "Test_ADQuadrature2D.h"

using namespace __Radiation;

/**
 * Evaluate the emission function.
 */
slibreal_t Test_ADQuadrature2D_Emission::Evaluate(
    RadiationParticle *rp, Vector<3> &n, slibreal_t sinMu, slibreal_t cosMu, bool pol
) {
    if (this->nwavelengths == 0)
        CalculateAngularDistribution(n, sinMu, cosMu);
    else {
        if (pol)
            CalculatePolarization(rp, n, sinMu, cosMu);
        else
            CalculateSpectrum(n, sinMu, cosMu);
    }

    return this->power;
}

/**
 * Angular distribution.
 */
void Test_ADQuadrature2D_Emission::CalculateAngularDistribution(
    Vector<3> &__UNUSED__(n), slibreal_t __UNUSED__(sinMu), slibreal_t __UNUSED__(cosMu)
) { this->power = 1.0; }

/**
 * Angular and spectral distribution.
 */
void Test_ADQuadrature2D_Emission::CalculateSpectrum(
    Vector<3> &__UNUSED__(n), slibreal_t __UNUSED__(sinMu), slibreal_t __UNUSED__(cosMu)
) {
    this->power = 1.0;
    for (unsigned int i = 0; i < nwavelengths; i++) {
        this->I[i] = (slibreal_t)(i+1);
    }
}

/**
 * Angular and spectral distribution with polarization.
 */
void Test_ADQuadrature2D_Emission::CalculatePolarization(
    RadiationParticle *__UNUSED__(rp),
    Vector<3> &__UNUSED__(n), slibreal_t __UNUSED__(sinMu), slibreal_t __UNUSED__(cosMu)
) {
    this->power = 1.0;
    for (unsigned int i = 0; i < nwavelengths; i++) {
        this->I[i] = (slibreal_t)(i+1);
        this->Q[i] = (slibreal_t)(i+1);
        this->U[i] = (slibreal_t)(i+1);
        this->V[i] = (slibreal_t)(i+1);
    }
}

