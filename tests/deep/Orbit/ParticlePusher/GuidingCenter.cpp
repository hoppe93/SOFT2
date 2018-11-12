/**
 * Test the guiding-center equations of motion.
 *
 * The following tests are implemented:
 *   [X] - Verify (resonable) energy and magnetic moment
 *         conservation in orbit solution.
 *   [ ] - Compare against SOFTv1 orbits.
 */

#include <softlib/config.h>
#include "Orbit/GuidingCenterEquation.h"
#include "Orbit/ParticlePusher.h"
#include "Test_ParticlePusher.h"

/**
 * Run all tests for the guiding-center equations
 * of motion.
 */
bool Test_ParticlePusher::GuidingCenterEquation() {
    bool success = true;

    if (!GuidingCenterEquation_Conservation(false))
        success = false;
    else
        this->PrintOK("Guiding-center solution (w/o drifts) conserves momentum and magnetic moment.");

    if (!GuidingCenterEquation_Conservation(true))
        success = false;
    else
        this->PrintOK("Guiding-center solution (w/ drifts) conserves momentum and magnetic moment.");

    if (!GuidingCenterEquation_TableCompare(true))
        success = false;
    else
        this->PrintOK("Guiding-center-with-drifts solution verified against SOFTv1 solution.");

    if (!GuidingCenterEquation_TableCompare(false))
        success = false;
    else
        this->PrintOK("Guiding-center-without-drifts solution verified against SOFTv1 solution.");

    return success;
}

/**
 * Test conservation properties of guiding-center
 * equations of motion.
 *
 * drifts: If true, solves the guiding-center equations
 *         of motion with drift terms.
 */
bool Test_ParticlePusher::GuidingCenterEquation_Conservation(bool drifts) {
    unsigned int i, ntau;
    slibreal_t mudiff, pdiff, *sol0, *sole;

    for (i = 0; i < ngcorbits; i++) {
        struct orbit *orb = guiding_center_orbits+i;

        Orbit *o = GenerateOrbit(orb->r, orb->p, orb->thetap, drifts);

        ntau = o->GetNTau();
        pdiff = fabs(o->GetP2(0) - o->GetP2(ntau-1)) / o->GetP2(0);
        sol0 = o->GetInternalSolution(0);
        sole = o->GetInternalSolution(ntau-1);

        mudiff = fabs(sol0[GuidingCenterEquation::COORD_MU]-sole[GuidingCenterEquation::COORD_MU]);

        if (pdiff > GC_TOLP) {
            this->PrintError("Guiding-center equation (%s drifts): momentum deviates by %e.", (drifts?"w/":"w/o"), pdiff);
            return false;
        }

        if (mudiff > GC_TOLMU) {
            this->PrintError("Guiding-center equation (%s drifts): magnetic moment deviates by %e.", (drifts?"w/":"w/o"), mudiff);
            return false;
        }
    }

    return true;
}

/**
 * Compare orbits to the tabulated SOFTv1 orbits.
 *
 * drifts: If true, solves the first order guiding-center equations
 *         of motion (with drift terms). Otherwise, solves the
 *         zeroth order guiding-center equations.
 */
bool Test_ParticlePusher::GuidingCenterEquation_TableCompare(bool drifts) {
    unsigned int i, j, nmax;
    slibreal_t *sol, dx, dy, dz, x, y, z;
    struct orbit *arr;

    if (drifts) {
        arr  = guiding_center_orbits;
        nmax = guiding_center_orbits_nmax;
    } else {
        arr  = guiding_center_orbits_nd;
        nmax = guiding_center_orbits_nd_nmax;
    }
    sol = new slibreal_t[6*nmax];

    for (i = 0; i < ngcorbits; i++) {
        struct orbit *orb = arr+i;

        ParticlePusher *pp = GenerateOrbit_push(orb->r, orb->p, orb->thetap, drifts);
        Integrator<6> *intg = pp->GetIntegrator1();

        intg->OutputDense(orb->elems, orb->t, sol);

        delete pp;

        for (j = 1; j < orb->elems; j++) {
            x = sol[j*6+0];
            y = sol[j*6+1];
            z = sol[j*6+2];

            if (x != 0) dx = fabs((orb->x[j]-x)/Rm);
            else dx = fabs(orb->x[j]-x);

            if (y != 0) dy = fabs((orb->y[j]-y)/Rm);
            else dy = fabs(orb->y[j]-y);

            if (z != 0) dz = fabs((orb->z[j]-z)/Rm);
            else dz = fabs(orb->z[j]-z);

            if (dx > GC_SOL_TOL) {
                this->PrintError("Guiding-center solution: X-component is incorrect (Delta: %e).",dx);
                return false;
            } else if (dy > GC_SOL_TOL) {
                this->PrintError("Guiding-center solution: Y-component is incorrect (Delta: %e).",dy);
                return false;
            } else if (dz > GC_SOL_TOL) {
                this->PrintError("Guiding-center solution: Z-component is incorrect. (Delta: %e)",dz);
                return false;
            }
        }
    }

    return true;
}
