/**
 * Wrapper of the dilogarithm.
 */

#include <gsl/gsl_sf_dilog.h>
#include <softlib/config.h>
#include "Tools/Radiation/Models/Cone/ConeBremsstrahlungEmission.h"

slibreal_t dilog_func(const slibreal_t x) {
    return gsl_sf_dilog(x);
}

