/**
 * Reverse projection.
 * Based on a model suggested by Ola Embréus and
 * developed by Mathias Hoppe.
 */

#include <algorithm>
#include <cmath>
#include <vector>
#include <softlib/Vector.h>
#include "Tools/Radiation/Models/Cone.h"
#include "Tools/Radiation/Models/Cone/Projection/Reverse.h"
#include "Tools/Radiation/RadiationParticle.h"

using namespace std;
using namespace __Radiation;

/**
 * Computes the fraction of the guiding-center cone
 * that overlaps the detector for the given particle
 * emission state.
 *
 * rp: Object representing the particle's emitting state.
 */
slibreal_t ConeProjectionReverse::ComputeOverlappingFraction(RadiationParticle *rp) {
    Vector<3> &phat = rp->GetPHat(),
        &e1 = detector->GetEHat1(),
        &e2 = detector->GetEHat2();
    Vector<3> g1, g2, r0, r;
    slibreal_t tanThetap, rd, t[8];
    unsigned int nt = 0, t0Intersects = 0;
    bool t0Inside;
    
    if (rp->GetPpar() != 0)
        tanThetap = fabs(rp->GetPperp() / rp->GetPpar());
    else {
        // TODO Implement radiation from stationary GCs
        return 0.0;
    }
    
    SetupBasisVectors(phat, g1, g2);

    rd = detector->GetAperture();
    /****************************
     * FIND INTERSECTION POINTS *
     ****************************/
    // Base edge of detector
    r0 = detector->GetPosition() + 0.5*rd*(e1 + e2) - rp->GetPosition();

    r=-rd*e1; FindCircleIntersectionPoints(r0, r, phat, g1, g2, tanThetap, t0Intersects, t, nt);
    r=-rd*e2; FindCircleIntersectionPoints(r0, r, phat, g1, g2, tanThetap, t0Intersects, t, nt);

    r0 -= rd*(e1 + e2);

    r=+rd*e1; FindCircleIntersectionPoints(r0, r, phat, g1, g2, tanThetap, t0Intersects, t, nt);
    r=+rd*e2; FindCircleIntersectionPoints(r0, r, phat, g1, g2, tanThetap, t0Intersects, t, nt);

    /************************************
     * CLASSIFY HIT (FULL/PARTIAL/MISS) *
     ************************************/
    // Is (x, y) = (tanThetap, 0) inside the detector?
    t0Inside = (t0Intersects%2 == 1);
    if (nt == 0) {  // No intersections?
        if (t0Inside) return 1.0;   // Full hit!
        else return 0.0;    // Miss
    }
    
    /****************************
     * SORT INTERSECTION POINTS *
     ****************************/
    vector<slibreal_t> tvec(t, t+nt);
    sort(tvec.begin(), tvec.end());

    /************************
     * SUM UP CONTRIBUTIONS *
     ************************/
    slibreal_t I = 0, s =-1.0;
    for (vector<slibreal_t>::iterator it = tvec.begin(); it != tvec.end(); it++) {
        I += s * (*it);
        s = -s;
    }

    if (t0Inside) return 1.0 - I;
    else return I;
}

/**
 * Find all intersection points between the emission circle
 * and the detector boundary described by the vector 'r',
 * originating from 'r0'.
 *
 * r0:           'Base' detector edge point
 * r:            Vector along which intersection should be
 *               checked (parametrisation is 'r0 + s*r',
 *               s in [0,1]).
 * phat:         Unit vector in direction of GC motion.
 * g1, g2:       Vectors spanning 'emission' plane
 *               (plane perpendicular to phat).
 * tanThetap:    Tangens of pitch angle.
 * t0Intersects: Number intersections between straight line
 *               along x from t0: (x, y) = (tanThetap, 0) and
 *               the detector edge under consideration.
 *               (variable updated by this function).
 * t:            Intersection points (updated by this function).
 * nt:           Number populated elements in 't' (updated
 *               by this function).
 */
void ConeProjectionReverse::FindCircleIntersectionPoints(
    Vector<3> &r0, Vector<3> &r, Vector<3> &phat,
    Vector<3> &g1, Vector<3> &g2,
    slibreal_t tanThetap, unsigned int &t0Intersects,
    slibreal_t t[8], unsigned int &nt
) {
    slibreal_t s, x, y, csd, tt;

    slibreal_t
        a1 = r0.Dot(g1),   a2 = r0.Dot(g2),
        b1 = r.Dot(g1),    b2 = r.Dot(g2),
        c  = r0.Dot(phat), d =  r.Dot(phat),
        tan2Thetap = tanThetap*tanThetap;

    slibreal_t
        den = b1*b1 + b2*b2 - d*d*tan2Thetap,
        num1 = c*d*tan2Thetap - a1*b1 - a2*b2,
        num2 = c*c*tan2Thetap - a1*a1 - a2*a2,
        num1den = num1 / den,
        sqr = num1den*num1den + num2 / den,
        rt;

    // Only allow GC to emit in the forward direction
    if (c+d <= 0)
        return;

    // Determine if the point
    //   (x,y) = (tanThetap, 0)
    // lies inside the detector
    if (b2 != 0) {
        s  = -a2/b2;
        tt = (a1*b2 - b1*a2) / (c*b2 - d*a2);

        if (0 <= s && s <= 1 &&
            tt >= tanThetap)
            t0Intersects++;
    }

    if (sqr < 0) return;
    else rt = sqrt(sqr);

    // Evaluate +- solutions
    if (num1den <= 1) { // + solution
        s = num1den + rt;
        
        csd = 1.0 / (c + s*d);
        x = (a1 + s*b1) * csd;
        y = (a2 + s*b2) * csd;
        if (0 <= s && s <= 1)
            t[nt++] = GetParameterT(x, y, tanThetap);
    }
    if (num1den >= 0) { // - solution
        s = num1den - rt;

        csd = 1.0 / (c + s*d);
        x = (a1 + s*b1) * csd;
        y = (a2 + s*b2) / (c + s*d);
        if (0 <= s && s <= 1)
            t[nt++] = GetParameterT(x, y, tanThetap);
    }
}

/**
 * Convert a cartesian point (x,y) to the corresponding
 * value of the parameter t used to parametrize the
 * emission circle.
 *
 * x:  Cartesian X position of the intersection point
 * y:  Cartesian Y position of the intersection point
 */
slibreal_t ConeProjectionReverse::GetParameterT(
    const slibreal_t x, const slibreal_t y, const slibreal_t tanThetap
) {
    slibreal_t p = x / tanThetap, t;
    if (p >= 1.0)
        t = 0.0;
    else if (p <= -1.0)
        t = 0.5;
    else
        t = acos(x / tanThetap) / (2.0*M_PI);

    if (y >= 0) return t;
    else return (1.0 - t);
}

/**
 * Setup the emission plane basis vectors g1 & g2.
 *
 * phat:   Unit vector in the direction of guiding-center motion.
 * g1, g2: On return, contains the new basis vectors spanning the emission plane.
 */
void ConeProjectionReverse::SetupBasisVectors(Vector<3> &phat, Vector<3> &g1, Vector<3> &g2) {
    if (phat[1] == 0) {
        g1[0] = g1[2] = 0;
        g1[1] = 1;
    } else {
        g1[0] = phat[1];
        g1[1] =-phat[0];

        g1.Normalize();
    }

    g2[0] = phat[1]*g1[2] - phat[2]*g1[1];
    g2[1] = phat[2]*g1[0] - phat[0]*g1[2];
    g2[2] = phat[0]*g1[1] - phat[1]*g1[0];
}

