/**
 * Helper functions for locating the SOV
 * in the toroidal integral.
 */

#include <functional>
#include <cmath>
#include "Tools/Radiation/sov_helpers.h"

#define SIGN(a,b)       ((b)<0?-fabs(a):fabs(a))

using namespace std;

namespace __Radiation {
/**
 * Shifts the given variables to the left so that
 *   a = b, b = c, c = d
 * on return.
 */
void shift3(slibreal_t &a, slibreal_t &b, slibreal_t &c, const slibreal_t d) {
    a = b;
    b = c;
    c = d;
}

/**
 * Computes "x mod 2*pi", i.e. ensures
 * that 'x' lies in the interval [0,2*pi].
 *
 * x: Value to compute modulo of.
 */
slibreal_t mod2pi(const slibreal_t x) {
    if (x < 0)
        return (x + 2.0*M_PI*ceil(-x/(2.0*M_PI)));
    else if (x > 2.0*M_PI)
        return (x - 2.0*M_PI*floor(x/(2.0*M_PI)));
    else
        return x;
}

/**
 * Bracket a minimum of the given function.
 * By setting 'ffac = -1' the maximum of the given
 * function can be bracketed instead.
 *
 * f:    Function to bracket optimum for.
 * ffac: Factor by which to multiply the given function.
 *       If +1, this function brackets the minimum, and
 *       if -1, it brackets the maximum.
 * ax:   On return, contains the lower limit of the
 *       interval in which the optimum must lie.
 * bx:   On return, contains the upper limit of the
 *       interval in which the optimum must lie.
 */
void bracket_min(function<slibreal_t(const slibreal_t)> f, const slibreal_t ffac, slibreal_t *axr, slibreal_t *bxr) {
    const slibreal_t dx = 0.1*M_PI;
    slibreal_t a = 0, b = dx, c = 2.0*dx;
    slibreal_t fa = ffac*f(a), fb = ffac*f(b), fc = ffac*f(c);

    while ((fa < fb || fb > fc) && a < 2.0*M_PI) {
        shift3(a, b, c, c+dx);
        shift3(fa, fb, fc, ffac*f(c));
    }

    *axr = a;
    *bxr = c;
}

/**
 * Brackets one root of the given function by
 * dividing the given interval into 'n' pieces and
 * checks if the root lies within each '2n'
 * subinterval.
 *
 * f:      Function to bracket roots for.
 * low:    Lower limit of interval to search.
 * high:   Upper limit of interval to search.
 * n:      Number of subintervals to partition
 *         the given interval into.
 * brackl: On return, contains the lower limit of the
 *         bracketed interval within which the root must lie.
 * brackr: On return, contains the upper limit of the
 *         bracketed interval within which the root must lie.
 *
 * RETURNS 'true' if a root was found in the given
 * interval. Otherwise returns 'false'.
 */
bool bracket_root(function<slibreal_t(const slibreal_t)> f, const slibreal_t low, const slibreal_t high, const unsigned int n, slibreal_t *brackl, slibreal_t *brackr) {
    slibreal_t dx = (high - low) / n;
    slibreal_t x = low, fp = f(x);

    for (unsigned int i = 0; i < n; i++) {
        slibreal_t fc = f(x += dx);
        if (fc*fp <= 0.0) {
            *brackl = x - dx;
            *brackr = x;
            return true;
        }

        fp = fc;
    }

    return false;
}

/**
 * Find the minimum of the given function.
 * By setting 'ffac = -1' the maximum of the given
 * function can be determined instead.
 *
 * f:    Function to optimize.
 * ffac: Factor to multiply function with. If '+1',
 *       the minimum of the function is located,
 *       and if '-1' then f's maximum is located.
 * tol:  Absolute tolerance within which to locate
 *       the optimum.
 */
slibreal_t find_min(function<slibreal_t(const slibreal_t)> f, const slibreal_t ffac, const slibreal_t tol) {
    slibreal_t a, b, c, d, yc, yd;

    bracket_min(f, ffac, &a, &b);

    slibreal_t
        h = b-a,
        invphi = 0.5*(sqrt(5.0) - 1.0),
        invphi2 = invphi*invphi;

    c = a + invphi2 * h;
    d = a + invphi  * h;
    yc = ffac*f(c);
    yd = ffac*f(d);

    while (h > tol) {
        if (yc < yd) {
            b = d;
            d = c;
            yd = yc;
            h = invphi * h;
            c = a + invphi2*h;
            yc = ffac*f(c);
        } else {
            a = c;
            c = d;
            yc = yd;
            h = invphi*h;
            d = a + invphi * h;
            yd = ffac*f(d);
        }
    }

    slibreal_t sol;
    if (yc < yd)
        sol = 0.5*(a+d);
    else
        sol = 0.5*(c+b);

    return mod2pi(sol);
}

/**
 * Use bisection to locate the root in the given
 * bracketed interval.
 *
 * f:   Function to locate root of.
 * a:   Lower limit of bracketed interval.
 * b:   Upper limit of bracketed interval.
 * tol: Absolute tolerance within which to
 *      locate the root.
 */
slibreal_t find_root(
    function<slibreal_t(const slibreal_t)> f,
    slibreal_t a, slibreal_t b, const slibreal_t tol
) {
    slibreal_t c, fa, fc;

    c = 0.5*(a+b);
    fa = f(a);
    fc = f(c);

    do {
        if (fa*fc >= 0) {
            a = c;
            fa = fc;
        } else {
            b = c;
        }

        c = 0.5*(a+b);
        fc = f(c);
    } while (fc != 0 && fabs(a-b) > tol);

    return mod2pi(0.5*(a+b));
}
};//namespace __Radiation
