
#include <gsl/gsl_machine.h>
#include <gsl/gsl_sf.h>
#include <softlib/config.h>

/*
 * Evaluates the ordinary Bessel function of the first kind
 */
slibreal_t cyclotron_func1(const unsigned int m, const slibreal_t x) {
	return (slibreal_t)(gsl_sf_bessel_Jn(m, x));
}

