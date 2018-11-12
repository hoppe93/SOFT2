/**
 * Implementation of the edge check condition.
 */

#include <softlib/Vector.h>
#include "Tools/Radiation/Models/Cone.h"
#include "Tools/Radiation/RadiationParticle.h"

using namespace __Radiation;

bool Cone::EdgeCheck(RadiationParticle *rp, const slibreal_t sinphi, const slibreal_t cosphi) {
    Vector<3> kappa = Vector<3>::Cross(rp->GetCurlBHat(), rp->GetBhat());

    slibreal_t kx = kappa[0], ky = kappa[1];

    kappa[0] = kx*cosphi + ky*sinphi;
    kappa[1] =-kx*sinphi + ky*cosphi;

    slibreal_t kRcp = kappa.Dot(rp->GetRCP());
    slibreal_t kE1  = kappa.Dot(this->parent->detector->GetEHat1());
    slibreal_t kE2  = kappa.Dot(this->parent->detector->GetEHat2());
    slibreal_t rd2  = 0.5 * this->parent->detector->GetAperture();

    /*slibreal_t
        d1 = kRcp - (rd2*kE1 + rd2*kE2),
        d2 = kRcp - (rd2*kE1 - rd2*kE2),
        d3 = kRcp + (rd2*kE1 + rd2*kE2),
        d4 = kRcp + (rd2*kE1 - rd2*kE2);*/
    slibreal_t
        r1 = rd2*fabs(kE1 + kE2),
        r2 = rd2*fabs(kE1 - kE2),
        rr = fabs(kRcp);

    return (rr < r1 || rr < r2);
}

