/**
 * Implementation of the emission class for the
 * ADQuadrature2D test module.
 */

#include <softlib/config.h>
#include <softlib/Vector.h>
#include "Test_ADQuadrature2D.h"

/**
 * Evaluate the emission function.
 */
slibreal_t Test_ADQuadrature2D_Emission::Evaluate(
    Vector<3> &n, slibreal_t sinMu, slibreal_t cosMu, bool pol
) {
    if (this->nwavelengths == 0)
        return CalculateAngularDistribution(n, sinMu, cosMu);
    else {
        if (pol)
            return CalculatePolarization(n, sinMu, cosMu);
        else
            return CalculateSpectrum(n, sinMu, cosMu);
    }
}

/**
 * Angular distribution.
 */
slibreal_t Test_ADQuadrature2D_Emission::CalculateAngularDistribution(
    Vector<3> &__UNUSED__(n), slibreal_t __UNUSED__(sinMu), slibreal_t __UNUSED__(cosMu)
) {
    this->power = 1.0;
    return this->power;
}

/**
 * Angular and spectral distribution.
 */
slibreal_t Test_ADQuadrature2D_Emission::CalculateSpectrum(
    Vector<3> &__UNUSED__(n), slibreal_t __UNUSED__(sinMu), slibreal_t __UNUSED__(cosMu)
) {
    this->power = 1.0;
    for (unsigned int i = 0; i < nwavelengths; i++) {
        this->I[i] = (slibreal_t)(i+1);
    }

    return this->power;
}

/**
 * Angular and spectral distribution with polarization.
 */
slibreal_t Test_ADQuadrature2D_Emission::CalculatePolarization(
    Vector<3> &__UNUSED__(n), slibreal_t __UNUSED__(sinMu), slibreal_t __UNUSED__(cosMu)
) {
    this->power = 1.0;
    for (unsigned int i = 0; i < nwavelengths; i++) {
        this->I[i] = (slibreal_t)(i+1);
        this->Q[i] = (slibreal_t)(i+1);
        this->U[i] = (slibreal_t)(i+1);
        this->V[i] = (slibreal_t)(i+1);
    }

    return this->power;
}

