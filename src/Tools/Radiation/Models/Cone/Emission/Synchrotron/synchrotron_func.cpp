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
 * (Based on [Aharonian et al. (2010; Appendix D)],
 *  https://doi.org/10.1103/PhysRevD.82.043002)
 */
/*slibreal_t synchrotron_func(const slibreal_t x) {
    slibreal_t x13 = cbrt(x),
               x23 = x13*x13,
               x43 = x23*x23;

    return 2.15*x13* sqrt(cbrt(1.0 + 3.06*x)) *
        (1.0 + 0.884*x23 + 0.471*x43) / (1.0 + 1.64*x23 + 0.974*x43) * exp(-x);
}
*/

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

