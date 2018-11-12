/**
 * Calculate the fraction of the cone that overlaps the detector.
 */

#include "Tools/Radiation/Models/Cone.h"
#include "Tools/Radiation/Models/Cone/Projection/Original.h"
#include "Tools/Radiation/RadiationParticle.h"

using namespace __Radiation;

/**
 * Computes the fraction of the guiding-center cone
 * that overlaps the detector for the given particle
 * emission state.
 *
 * rp: Object representing the particle's emitting state.
 */
slibreal_t ConeProjectionOriginal::ComputeOverlappingFraction(RadiationParticle *rp) {
    slibreal_t a, b, x0, X, cosThetap, sinThetap, cosPhi, sinPhi, p,
        XcosPhi, ccms, ppar, pperp, sinXi=0, cosXi=1, xh, yh, xoff, yoff,
        rd_2, forwardCone;

    Vector<3> &e1 = detector->GetEHat1(),
              &e2 = detector->GetEHat2(),
              &nhat = detector->GetDirection(),
              &phat = rp->GetPHat(),
              &rcp = rp->GetRCP(),
              xhit;

    ppar = rp->GetPpar();
    pperp = rp->GetPperp();
    p = sqrt(rp->GetP2());

    cosThetap = fabs(ppar / p);
    sinThetap = pperp / p;
    cosPhi = phat.Dot(nhat);
    sinPhi = sqrt(1 - cosPhi*cosPhi);
    X = fabs(rcp.Dot(nhat) / cosPhi);
    rd_2 = 0.5 * detector->GetAperture();

    if (sinPhi != 0) {
        sinXi = phat.Dot(e2) / sinPhi;
        cosXi = phat.Dot(e1) / sinPhi;
    }

    ccms = cosThetap*cosThetap - sinPhi*sinPhi;
    XcosPhi = X*cosPhi;

    a  = fabs(cosThetap*sinThetap / ccms * XcosPhi);
    b  = fabs(sinThetap / sqrt(fabs(ccms)) * XcosPhi);
    x0 = -sinThetap*sinThetap / ccms * X*sinPhi;

    xhit = X*phat + rcp;
    xh = xhit[0]*nhat[1] + xhit[1]*e1[1] + xhit[2]*e2[1];
    yh = xhit[0]*nhat[2] + xhit[1]*e1[2] + xhit[2]*e2[2];

    xoff = -x0*cosXi + xh;
    yoff =  x0*sinXi + yh;

    forwardCone = -sinThetap*sinPhi / ccms * (cosThetap*cosPhi + sinThetap*sinPhi);
    if (cosPhi < 0) {
        if (forwardCone > 1.0) return 0;
    } else if (forwardCone < 1.0 || ccms > 0) return 0;

    struct intersection_points ip;
    if (ccms > 0) { // Ellipse
        EllipseGetIntersections(a, b, x0, sinXi, cosXi, xh, yh, &ip);
    } else if (ccms == 0) { // Parabola
        // not handled...
        return 0.0;
    } else {    // Hyperbola
        HyperbolaGetIntersections(a, b, x0, sinXi, cosXi, xh, yh, &ip);
    }

    if (ip.nt < 2) {
        if (ccms > 0) {     /* If ellipse... */
            slibreal_t x = xoff + a*cosXi;
            slibreal_t y = yoff - a*sinXi;

            if (fabs(x) < rd_2 && fabs(y) < rd_2) {    /* Full hit? */
                return 1.0;
            } else return 0.0;      /* No hit */
        } else return 0.0;  /* Hyperbola's cannot be fully contained within detector */
    }

    /***************************************
     * COMPUTE CORRESPONDING POINT ON CONE *
     ***************************************/
    SortIntersectionPoints(ip);

    slibreal_t sgn = -1.0,
        pnt = 0.5*(ip.t[1]+ip.t[0]),
        lastpoint = 0.0,
        fraction = 0.0, contrib,
        xp, yp;
    int count = ip.nt;

    // Prepare for the transformation
    InitReverseProjection(rp->GetPHat(), sinXi, cosXi);

    /* QUICK REMINDER
     * 'From'-points are points you INTEGRATE from.
     * 'To'-points are points you INTEGRATE to.
     */

    // If an ellipse...
    if (ccms > 0) {
        // If the first point is a 'to'-point, change sign
        if (fabs(xoff+a*cosXi*cos(pnt)+b*sinXi*sin(pnt)) >= rd_2 ||
            fabs(yoff-a*sinXi*cos(pnt)+b*cosXi*sin(pnt)) >= rd_2) {
            sgn = 1;

            // Get its corresponding 'from'-point.
            xp = ip.x[ip.nt-1] - xoff,
            yp = ip.y[ip.nt-1] - yoff;

            lastpoint = ToPointOnUnitCone(xp, yp);
            // No need to handle the final point again
            count--;
            fraction += -lastpoint;
        }
    } /* Note: the first pint of a hyperbola MUST be a 'from'-point,
         since its domain is (-inf, inf). */

    for (int i = 0; i < count; i++) {
        xp = ip.x[i] - xoff;
        yp = ip.y[i] - yoff;

        contrib = ToPointOnUnitCone(xp, yp);

        fraction += sgn*contrib;
        /* If this is a 'to'-point and this point
         * is smaller than the previous pint, the
         * line segment connecting them passes
         * through 0. Therefore, add 2*PI. */
        if (sgn > 0 && contrib < lastpoint)
            fraction += 2.0*M_PI;

        lastpoint = contrib;
        sgn = -sgn;
    }

    fraction /= 2.0*M_PI;
    return fraction;
}

/**
 * Bubble sort the intersection points in order
 * from smallest to largest (array will be like
 * [0 -> 2*PI] or [-inf -> inf]).
 *
 * ip: Intersection point container to sort.
 */
void ConeProjectionOriginal::SortIntersectionPoints(struct intersection_points &ip) {
     int i, swapped=1, tlen=ip.nt;
     slibreal_t swp;

     while (swapped) {
        swapped = 0;
        for (i = 1; i < tlen; i++) {
            if (ip.t[i-1] > ip.t[i]) {
                swp = ip.t[i];
                ip.t[i] = ip.t[i-1];
                ip.t[i-1] = swp;

                swp = ip.x[i];
                ip.x[i] = ip.x[i-1];
                ip.x[i-1] = swp;

                swp = ip.y[i];
                ip.y[i] = ip.y[i-1];
                ip.y[i-1] = swp;

                swapped = 1;
            }
        }
        tlen--;
     }
}

/**
 * Prepare for transforming from the detector plane
 * to the unit cone, in order to compute the total
 * overlapping fraction of the guiding-center cone.
 *
 * vhat:  Guiding-center direction of motion (V / |V|).
 * sinxi: Sine of rotation angle of projected curve
 *        in the detector plane.
 * cosxi: Cosine of rotation angle of projected curve
 *        in the detector plane.
 */
void ConeProjectionOriginal::InitReverseProjection(
    Vector<3> &vhat, slibreal_t sinxi, slibreal_t cosxi
) {
    Vector<3> &e1 = detector->GetEHat1(),
              &e2 = detector->GetEHat2();
    int i, j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            if (i == j) Pn[i][j] = 1.0;
            else Pn[i][j] = 0.0;

            Pn[i][j] -= vhat[i]*vhat[j];
        }
    }

    e1r = e1*cosxi + e2*sinxi;
    e2r =-e1*sinxi + e2*cosxi;

    Pne1r[0] = e1r.Dot(Pn[0]);
    Pne1r[1] = e1r.Dot(Pn[1]);
    Pne1r[2] = e1r.Dot(Pn[2]);

    Pne2r[0] = e2r.Dot(Pn[0]);
    Pne2r[1] = e2r.Dot(Pn[1]);
    Pne2r[2] = e2r.Dot(Pn[2]);

    absPne1r = Pne1r.Norm();
    absPne2r = Pne2r.Norm();
}

/**
 * Convert a point along the projected curve to
 * a point on the unit cone.
 *
 * x, y: X & Y coordinates of the intersection point
 *       (in the system which has inclination along
 *       the x-axis only).
 */
slibreal_t ConeProjectionOriginal::ToPointOnUnitCone(
    slibreal_t x, slibreal_t y
) {
    Vector<3> r, Pnx,
        &e1 = detector->GetEHat1(),
        &e2 = detector->GetEHat2();

    r = x*e1 + y*e2;

    Pnx[0] = r.Dot(Pn[0]);
    Pnx[1] = r.Dot(Pn[1]);
    Pnx[2] = r.Dot(Pn[2]);

    slibreal_t
        absPnx = Pnx.Norm(),
        cosDelta = Pnx.Dot(Pne1r) / (absPnx*absPne1r),
        sinDelta = Pnx.Dot(Pne2r) / (absPnx*absPne2r);

    if (sinDelta > 0)
        return acos(cosDelta);
    else
        return 2.0*M_PI - acos(cosDelta);
}

/**************************************
 * ELLIPSE                            *
 **************************************/
/**
 * Get the points on an ellipse that intersect
 * the detector.
 *
 * a, b, x0: Ellipse parameters, as defined in the
 *           SOFT manual.
 * ip:       Contains information about the intersection
 *           points upon return.
 */
void ConeProjectionOriginal::EllipseGetIntersections(
    slibreal_t a, slibreal_t b, slibreal_t x0,
    slibreal_t sinXi, slibreal_t cosXi,
    slibreal_t xhit, slibreal_t yhit,
    struct intersection_points *ip
) {
    slibreal_t Ax, Ay, Bx, By, dx, dy, rd_2;

    Ax = a*cosXi,  Ay = a*sinXi;
    Bx = b*sinXi,  By = b*cosXi;
    dx = -x0*cosXi + xhit;
    dy =  x0*sinXi + yhit;
    rd_2 = 0.5 * detector->GetAperture();

    ip->nt = 0;
    // x = +rd/2,  |y| < rd/2
    EllipseGetIntersection(Ax, Ay, Bx, By, +rd_2, dx, dy, ip->x, ip->y, ip);

    // x = -rd/2,  |y| < rd/2
    EllipseGetIntersection(Ax, Ay, Bx, By, -rd_2, dx, dy, ip->x, ip->y, ip);

    // y = +rd/2,  |x| < rd/2
    EllipseGetIntersection(-Ay, -Ax, By, Bx, +rd_2, dy, dx, ip->y, ip->x, ip);

    // y = -rd/2,  |x| < rd/2
    EllipseGetIntersection(-Ay, -Ax, By, Bx, -rd_2, dy, dx, ip->y, ip->x, ip);
}

/**
 * Find intersection points (internal). Evaluates (6.17)
 * in the SOFT manual.
 * 
 * A, B, d: Coefficients of parametrisation (see manual).
 * rd_2:    Half detector length (with sign, depending on which
 *          detector corner is being considered).
 * flipped: X and Y coordinates should be flipped in ip.
 * ip:      Contains intersection points on return.
 */
void ConeProjectionOriginal::EllipseGetIntersection(
    slibreal_t Ax, slibreal_t Ay,
    slibreal_t Bx, slibreal_t By,
    slibreal_t rd_2,
    slibreal_t d, slibreal_t origin,
    slibreal_t x[8], slibreal_t y[8],
    struct intersection_points *ip
) {
    slibreal_t d2 = d - rd_2;
    slibreal_t dif = Ax - d2;       // XXX why the minus?
    if (dif == 0.0) return;
    
    slibreal_t fac = Ax*Ax + Bx*Bx - d2*d2;
    if (fac < 0.0) return;

    slibreal_t a = sqrt(fac) / dif,
               b = Bx / dif;
    slibreal_t fac1 = b+a, fac12 = fac1*fac1;
    slibreal_t fac2 = b-a, fac22 = fac2*fac2;
    slibreal_t sqr1 = 1.0 / (1.0 + fac12);
    slibreal_t sqr2 = 1.0 / (1.0 + fac22);

    slibreal_t t1 = 2.0*atan(fac1);
    slibreal_t t2 = 2.0*atan(fac2);

    slibreal_t cost1 = (1.0 - fac12) * sqr1,
               cost2 = (1.0 - fac22) * sqr2,
               sint1 = 2.0*fac1*sqr1,
               sint2 = 2.0*fac2*sqr2;
    
    slibreal_t y1 = origin - Ay*cost1 + By*sint1;
    if (fabs(y1) <= fabs(rd_2)) {
        ip->t[ip->nt] = t1;
        x[ip->nt] = rd_2;
        y[ip->nt] = y1;
        ip->nt++;
    }

    slibreal_t y2 = origin - Ay*cost2 + By*sint2;
    if (fabs(y2) <= fabs(rd_2)) {
        ip->t[ip->nt] = t2;
        x[ip->nt] = rd_2;
        y[ip->nt] = y2;
        ip->nt++;
    }
}

/**************************************
 * HYPERBOLA                          *
 **************************************/
/**
 * Get the points on a hyperbola that intersect
 * the detector.
 *
 * a, b, x0: Hyperbola parameters, as defined in the
 *           SOFT manual.
 * ip:       Contains information about the intersection
 *           points upon return.
 */
void ConeProjectionOriginal::HyperbolaGetIntersections(
    slibreal_t a, slibreal_t b, slibreal_t x0,
    slibreal_t sinXi, slibreal_t cosXi,
    slibreal_t xhit, slibreal_t yhit,
    struct intersection_points *ip
) {
    slibreal_t Ax, Ay, Bx, By, dx, dy, rd_2;

    Ax = a*cosXi,  Ay = a*sinXi;
    Bx = b*sinXi,  By = b*cosXi;
    dx = -x0*cosXi + xhit;
    dy =  x0*sinXi + yhit;
    rd_2 = 0.5 * detector->GetAperture();

    ip->nt = 0;
    // x = +rd/2,  |y| < rd/2
    HyperbolaGetIntersection(Ax, Ay, Bx, By, +rd_2, dx, dy, ip->x, ip->y, ip);

    // x = -rd/2,  |y| < rd/2
    HyperbolaGetIntersection(Ax, Ay, Bx, By, -rd_2, dx, dy, ip->x, ip->y, ip);

    // y = +rd/2,  |x| < rd/2
    HyperbolaGetIntersection(-Ay, -Ax, By, Bx, +rd_2, dy, dx, ip->y, ip->x, ip);

    // y = -rd/2,  |x| < rd/2
    HyperbolaGetIntersection(-Ay, -Ax, By, Bx, -rd_2, dy, dx, ip->y, ip->x, ip);
}

/**
 * Find intersection points (internal). Evaluates (6.17)
 * in the SOFT manual.
 * 
 * A, B, d: Coefficients of parametrisation (see manual).
 * rd_2:    Half detector length (with sign, depending on which
 *          detector corner is being considered).
 * flipped: X and Y coordinates should be flipped in ip.
 * ip:      Contains intersection points on return.
 */
void ConeProjectionOriginal::HyperbolaGetIntersection(
    slibreal_t Ax, slibreal_t Ay,
    slibreal_t Bx, slibreal_t By,
    slibreal_t rd_2,
    slibreal_t d, slibreal_t origin,
    slibreal_t x[8], slibreal_t y[8],
    struct intersection_points *ip
) {
    slibreal_t u, b, sqr, arg1, arg2, d2, t, yval,
        sinht, cosht;
    if (Ax+Bx == 0) return;

    d2 = d-rd_2;

    u = d2*d2 - Ax*Ax + Bx*Bx;
    b = 1/(Ax+Bx);

    if (u < 0) return;

    sqr = sqrt(u);
    arg1 = (sqr-d2)*b;
    arg2 = (-sqr-d2)*b;

    if (arg1 > 0) {
        t = log(arg1);
        cosht = cosh(t), sinht = sinh(t);
        yval = origin - Ay*cosht + By*sinht;

        if (fabs(yval) <= fabs(rd_2)) {
            ip->t[ip->nt] = t;
            x[ip->nt] = rd_2;
            y[ip->nt] = yval;
            ip->nt++;
        }
    }

    if (arg2 > 0) {
        t = log(arg2);
        cosht = cosh(t), sinht = sinh(t);
        yval = origin - Ay*cosht + By*sinht;

        if (fabs(yval) <= fabs(rd_2)) {
            ip->t[ip->nt] = t;
            x[ip->nt] = rd_2;
            y[ip->nt] = yval;
            ip->nt++;
        }
    }
}

