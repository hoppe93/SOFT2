/**
 * Implementation of basic methods in SOFTEquation.
 */

#include <omp.h>
#include <softlib/config.h>
#include <softlib/Vector.h>
#include "Orbit/SOFTEquation.h"

/**
 * Get particle radial position.
 *
 * x: Solution vector (6-dimensional).
 */
slibreal_t SOFTEquation::GetPositionR(const slibreal_t *x) {
    return GetPositionR(
        x[0], x[1], x[2],
        x[3], x[4], x[5]
    );
}
slibreal_t SOFTEquation::GetPositionR(const Vector<6>& x) {
    return GetPositionR(
        x[0], x[1], x[2],
        x[3], x[4], x[5]
    );
}

/**
 * Get particle vertical position.
 * 
 * x: Solution vector (6-dimensional).
 */
slibreal_t SOFTEquation::GetPositionZ(const slibreal_t *x) {
    return GetPositionZ(
        x[0], x[1], x[2],
        x[3], x[4], x[5]
    );
}
slibreal_t SOFTEquation::GetPositionZ(const Vector<6>& x) {
    return GetPositionZ(
        x[0], x[1], x[2],
        x[3], x[4], x[5]
    );
}

/**
 * Calculates the phases-space Jacobians denoted
 * 'J' and 'Jp' in the SOFT documentation from
 * the other quantities of this SOFTEquation.
 *
 * solution:  Solution from first integrator.
 * solution2: Solution from second integrator, launched
 *            a distance 'nudge' from the first.
 * o:         Orbit object to store result in.
 * nudge:     Nudge value used when calculating 'solution2'.
 */
void SOFTEquation::CalculateJacobians(slibreal_t *solution, slibreal_t *solution2, Orbit *o, slibreal_t nudge) {
    unsigned int i, ti, nt = o->GetNTau();

    slibreal_t
        *p = o->GetP(),
        *Babs = o->GetBabs(),
        *ppar = o->GetPpar(),
        *gamma = o->GetGamma(),
        *Jdtdrho = o->GetJdtdrho(),
        *Jp = o->GetJp(),
        dtau = o->GetTau()[1]-o->GetTau()[0];

    for (i = 0; i < nt; i++) {
        // Calculate spatial Jacobian?
        if (solution2 != nullptr) {
            slibreal_t dR_dt, dZ_dt, dR_drho, dZ_drho,
                X1, Y1, Z1, X2, Y2, Z2, R1;

            ti = i*3;
            X1 = solution [i*6+0];
            Y1 = solution [i*6+1];
            Z1 = solution [i*6+2];
            X2 = solution2[i*6+0];
            Y2 = solution2[i*6+1];
            Z2 = solution2[i*6+2];
            R1 = hypot(X1, Y1);

            dR_dt = (X1*p[ti] + Y1*p[ti+1]) * LIGHTSPEED / (R1*gamma[i]) * dtau;
            dZ_dt = p[ti+2] * LIGHTSPEED / gamma[i] * dtau;

            dR_drho = hypot(X2, Y2) - R1;
            dZ_drho = Z2 - Z1;

            Jdtdrho[i] = fabs(dR_dt*dZ_drho - dZ_dt*dR_drho) / nudge;
            if (this->particle != 0)
                Jdtdrho[i] *= this->particle->GetDRho();
        }

        Jp[i] = Babs[i] / Babs[0] * fabs(ppar[0] / ppar[i]);
    }
}

