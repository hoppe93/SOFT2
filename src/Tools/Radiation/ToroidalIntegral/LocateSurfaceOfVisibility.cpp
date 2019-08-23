/**
 * Find intersection points.
 */

#include <functional>
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/sov_helpers.h"

#define MAX(x,y)        ((x)>(y)?(x):(y))
#define MIN(x,y)        ((x)<(y)?(x):(y))

using namespace __Radiation;
using namespace std;

/**
 * Locate the toroidal angles at which the surface
 * of visibility for the given particle can be seen.
 * This method solves the SOV equation numerically
 * for the one or two toroidal angles at which it
 * may be visible.
 *
 * rp:   RadiationParticle object describing the
 *       particle emitting state.
 * phi1: First toroidal angle at which the SOV
 *       can be observed.
 * phi2: Second toroidal angle at which the SOV
 *       can be observed.
 *
 * NOTE  This function always returns a toroidal
 *       angle in 'phi1' that should be tested.
 *       If 'phi1' and 'phi2' are the same, this
 *       indicates a double root and only one of
 *       them should be tested.
 */
void Radiation::LocateSurfaceOfVisibility(RadiationParticle *rp, unsigned int *phi1, unsigned int *phi2) {
    const slibreal_t tol = 0.5*M_PI / ntoroidal;    // (2*pi/ntoroidal) / 2
    Vector<3> v = rp->GetPHat(), x = rp->GetPosition(), d = this->detector->GetPosition();

    slibreal_t
        vx = v[0], vy = v[1], vz = v[2],
        x0 = x[0], y0 = x[1], z0 = x[2],
        xd = d[0], yd = d[1], zd = d[2],
        // Here we let 'v' contain the direction of the radiation,
        // and only use the magnitude of 'cosThetap'
        cosThetap = fabs(rp->GetPpar()) / sqrt(rp->GetP2());

    /**********************************************
     * FUNCTION FOR EVALUATING V. R - cos(thetap) *
     **********************************************/
    function<slibreal_t(const slibreal_t)> f = [&vx,&vy,&vz,&x0,&y0,&z0,&xd,&yd,&zd,&cosThetap](const slibreal_t t) {
        slibreal_t cost = cos(t), sint = sin(t);
        slibreal_t rx, ry, rz, vhx, vhy, vhz, r, vr;

        rx = xd - x0*cost - y0*sint;
        ry = yd + x0*sint - y0*cost;
        rz = zd - z0;
        r  = sqrt(rx*rx + ry*ry + rz*rz);

        vhx = vx*cost + vy*sint;
        vhy =-vx*sint + vy*cost;
        vhz = vz;

        vr = (rx*vhx + ry*vhy + rz*vhz) / r;
        return (vr-cosThetap);
    };

    // Bracket roots
    slibreal_t l, h, t1, t2;
    unsigned i, j;

    /******************
     * LOCATE MAXIMUM *
     ******************/
    slibreal_t mx = find_min(f, -1, tol), mn;

    /*****************
     * BRACKET ROOTS *
     *****************/
    // No root (or possibly double root)
    if (f(mx) < 0) {
        t1 = t2 = mx;
    } else {
        mn = find_min(f, +1, tol);
        l = MIN(mx, mn);
        h = MAX(mx, mn);

        t1 = find_root(f, l, h, tol);
        t2 = find_root(f, h-2.0*M_PI, l, tol);
    }

    i = (unsigned int)round((ntoroidal-1)*t1/(2.0*M_PI));
    j = (unsigned int)round((ntoroidal-1)*t2/(2.0*M_PI));

    *phi1 = (i <= j ? i : j);
    *phi2 = (i >  j ? i : j);
}

unsigned int Radiation::LocatePointOfVisibility(RadiationParticle *rp) {
    Vector<3> v = rp->GetPHat(), x = rp->GetPosition(), d = this->detector->GetPosition();
    slibreal_t
        vx = v[0], vy = v[1],
        xp = x[0], yp = x[1],
        x0 = d[0], y0 = d[1];

        
    slibreal_t par1 = vx*y0-vy*x0, 
        par2 = vy*y0 + vx*x0,
        par3 = vy*xp - vx*yp;
    slibreal_t A = par1*par1 + par2*par2,
        B = par1*par3,
        C = par3*par3 - par2*par2;

    slibreal_t cosphi = (-B + sqrt(B*B-A*C))/A,
        sinphi = sqrt(1-cosphi*cosphi),
        a = (y0 + xp*sinphi-yp*cosphi)/(vy*cosphi - vx*sinphi),
        tol = 1e-7,
        phi1,
        eq_1 = (xp+a*vx)*cosphi + (yp+a*vy)*sinphi - x0,
        eq_2 = -(xp+a*vx)*sinphi + (yp+a*vy)*cosphi - y0;

    //phi1 = acos(cosphi);
    //return (unsigned int)round((ntoroidal-1)*phi1/(2.0*M_PI));
    
    if (fabs(eq_1) < tol && fabs(eq_2) < tol){
        phi1 = acos(cosphi);
        return (unsigned int)round((ntoroidal-1)*phi1/(2.0*M_PI));
    }
    
    sinphi = -sinphi;
    a = (y0 + xp*sinphi-yp*cosphi)/(vy*cosphi - vx*sinphi);
    eq_1 = (xp+a*vx)*cosphi + (yp+a*vy)*sinphi - x0;
    eq_2 = -(xp+a*vx)*sinphi + (yp+a*vy)*cosphi - y0;
    //printf("Second, abs eq1 = %e, abs eq2 = %e \n", fabs(eq_1), fabs(eq_2));
    if (fabs(eq_1) < tol && fabs(eq_2) < tol){
        phi1 = acos(cosphi) + M_PI;
        return (unsigned int)round((ntoroidal-1)*phi1/(2.0*M_PI));
    }
    
    cosphi = (-B - sqrt(B*B-A*C))/A;
    sinphi = sqrt(1-cosphi*cosphi);
    a = (y0 + xp*sinphi-yp*cosphi)/(vy*cosphi - vx*sinphi);
    eq_1 = (xp+a*vx)*cosphi + (yp+a*vy)*sinphi - x0;
    eq_2 = -(xp+a*vx)*sinphi + (yp+a*vy)*cosphi - y0;
    //printf("Third, abs eq1 = %e, abs eq2 = %e \n", fabs(eq_1), fabs(eq_2));
    if (fabs(eq_1) < tol && fabs(eq_2) < tol){
        phi1 = acos(cosphi);
        return (unsigned int)round((ntoroidal-1)*phi1/(2.0*M_PI));
    } 
    
    sinphi = -sinphi;
    a = (y0 + xp*sinphi-yp*cosphi)/(vy*cosphi - vx*sinphi);
    eq_1 = (xp+a*vx)*cosphi + (yp+a*vy)*sinphi - x0;
    eq_2 = -(xp+a*vx)*sinphi + (yp+a*vy)*cosphi - y0;
    //printf("Last, abs eq1 = %e, abs eq2 = %e \n", fabs(eq_1), fabs(eq_2));
    if (fabs(eq_1) < tol && fabs(eq_2) < tol){
        phi1 = acos(cosphi) + M_PI;
        return (unsigned int)round((ntoroidal-1)*phi1/(2.0*M_PI));
    }    
    
    printf("No solution \n");
    return 0;

        

/*
    slibreal_t cosphi_p = (-B + sqrt(B*B-A*C))/A,
        cosphi_m = (-B - sqrt(B*B-A*C))/A,
        sinphi_p = sqrt(1-cosphi_p*cosphi_p),
        a = (y0 + xp*sinphi_p-yp*cosphi_p)/(vy*cosphi_p - vx*sinphi_p),
        phi1;
    
    if (a > 0)
        phi1 = acos(cosphi_p);
    else
        phi1 = acos(cosphi_m);

    return (unsigned int)round((ntoroidal-1)*phi1/(2.0*M_PI));
*/

}

