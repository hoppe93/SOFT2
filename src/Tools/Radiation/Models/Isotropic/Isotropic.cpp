/**
 * Implementation of the 'Isotropic' radiation model.
 */

#include <cstdio>
#include <cmath>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "Orbit/Orbit.h"
#include "Tools/Radiation/Models/Isotropic.h"

using namespace __Radiation;

/**
 * Calculate isotropic radiation from the given
 * particle.
 *
 * rp: Particle state to calculate radiation for.
 */
void Isotropic::HandleParticle(RadiationParticle *rp, orbit_type_t otype, const slibreal_t, const slibreal_t) {
    Vector<3> &rcp = rp->GetRCP(),
              &nhat = parent->detector->GetDirection(),
              &ehat1 = parent->detector->GetEHat1(),
              &ehat2 = parent->detector->GetEHat2();

    slibreal_t cn = nhat.Dot(rcp),
               c1 = ehat1.Dot(rcp),
               c2 = ehat2.Dot(rcp),
               aperture = 0.5*parent->detector->GetAperture();

    slibreal_t
        App = LambdaIso(+aperture, +aperture, cn, c1, c2),
        Amm = LambdaIso(-aperture, -aperture, cn, c1, c2),
        Amp = LambdaIso(-aperture, +aperture, cn, c1, c2),
        Apm = LambdaIso(+aperture, -aperture, cn, c1, c2);

    this->power = this->value * (App + Amm - Amp - Apm);
}

/**
 * The isotropic emission function (see the SOFT manual,
 * chapter "Isotropic radiation") \Lambda_{\mathrm{iso}}(l1,l2; x-X_0)
 */
slibreal_t Isotropic::LambdaIso(
    slibreal_t l1, slibreal_t l2,
    slibreal_t cn, slibreal_t c1, slibreal_t c2
) {
    slibreal_t arg =
        (c1-l1)*(c2-l2) / (cn*sqrt(cn*cn + (c1-l1)*(c1-l1) + (c2-l2)*(c2-l2)));

    return atan(arg);
}

