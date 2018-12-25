/**
 * Implementation of the first synchrotron function.
 */

#include <gsl/gsl_machine.h>
#include <gsl/gsl_sf.h>
#include <softlib/config.h>
#include "Tools/Radiation/Models/Cone/ConeSynchrotronEmission.h"

/**
 * Evaluates the first synchrotron function:
 *
 *   F(x) = x * Integral[ K_{5/3}(t), {t, x, inf}]
 *
 * (Based on the implementation in GSL)
 */
slibreal_t synchrotron_func1(const slibreal_t x) {
    if (x >= -8.0*GSL_LOG_DBL_MIN/7.0)
        return 0.0; /* Beyond machine precision */
    else
        return (slibreal_t)gsl_sf_synchrotron_1(x);
}

/**
 * Evaluates the second synchrotron function:
 *
 *   F(x) = x * K_{2/3}(x)
 *
 * (Based on the implementation in GSL)
 */
slibreal_t synchrotron_func2(const slibreal_t x) {
    if (x >= -8.0*GSL_LOG_DBL_MIN/7.0)
        return 0.0; /* Beyond machine precision */
    else
        return (slibreal_t)gsl_sf_synchrotron_2(x);
}

/**
 * Evaluates the "third" synchrotron function (our terminology):
 *
 *   F(x) = x * K_{1/3}(x)
 *
 * (Based on implementation in GSL)
 */
slibreal_t synchrotron_func3(const slibreal_t x) {
    return (slibreal_t)(x*gsl_sf_bessel_Knu(1.0/3.0, x));
}

